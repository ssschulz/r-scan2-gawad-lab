# Read in each SCAN2 object and retain sites that would pass if the
# calling threshold was increased to target.fdr=0.5 (which is very high).
# The entire object list can't be retained in memory since it requires
# 2-3 GB of RAM per sample.
#
# Future optimization of the table might make this feasible.
prepare.object <- function(object, quiet=FALSE) {
    newgatk <- object@gatk[static.filter & lysis.fdr <= 0.5 & mda.fdr <= 0.5]
    if (!quiet) {
        cat(sprintf("%s: considering %d candidates (%d high quality passing calls)\n",
            object@single.cell, nrow(newgatk), nrow(object@gatk[pass == TRUE])))
    }
    object@gatk <- newgatk  # XXX: likely need a formal way of replacing this table with a subset..
    compute.filter.reasons(object)
    object
}


get.objects.for.sig.rescue <- function(object.paths, quiet=FALSE) {
    os <- lapply(object.paths, function(path) {
        if (!quiet) print(path)
        x <- get(load(path))
        prepare.object(x)
    })
    names(os) <- sapply(os, function(o) o@single.cell)
    os
}


# Modifies 'o' by reference. The returned o is not a particularly valid
# SCAN2 object.
#
# N.B. there are currently no reasons to separate SNVs and indels here
# because the various test columns already incoporate differences in
# calling parameters.
compute.filter.reasons <- function(o, target.fdr=o@call.mutations$target.fdr) {
    m <- o@gatk[, .(pass=!pass, abc.test,
        lysis.test=!is.na(lysis.fdr) & lysis.fdr <= target.fdr,
        mda.test=!is.na(mda.fdr) & mda.fdr <= target.fdr,
        cigar.id.test, cigar.hs.test,
        lowmq.test, dp.test, min.sc.alt.test,
        csf.test, dbsnp.test)]

    # some tests can be NA - e.g., many tests when dp=0, cross-sample filter test when
    # the site is not in the panel, etc.  these tests should be considered to
    # have failed when NA.
    m[is.na(m)] <- FALSE

    # don't overlook the negation apply(!m, ...
    o@gatk[, filter.reasons :=
        apply(!m, 1, function(row) paste(colnames(m)[row], collapse='&'))]
}


setGeneric("mutsig.rescue.one", function(object, artifact.sig, true.sig, rescue.target.fdr=0.01, muttype=c('snv', 'indel'))
    standardGeneric("mutsig.rescue.one"))
setMethod("mutsig.rescue.one", "SCAN2", function(object, artifact.sig, true.sig, rescue.target.fdr=0.01, muttype=c('snv', 'indel')) {
    mt <- match.arg(muttype)

    # All work in this function will be done on a copy of the object with a much, much
    # smaller GATK table.  Results will be joined back at the end.
    tmpo <- prepare.object(object)

    sigtype <- if (mt == 'snv') sbs96 else id83
    mutsigs <- sigtype(tmpo@gatk[muttype == mt & filter.reasons == 'lysis.test']$mutsig)

    sigscores <- get.sig.score(mutsigs=mutsigs,
        artifact.sig=artifact.sig, true.sig=true.sig)

    # it doesn't seem to be possible to use a column assigned by := for another
    # assignment in the same data.table statement.  i.e., to combine all of these
    # into a single statement.
    tmpo@gatk[muttype == mt & filter.reasons == 'lysis.test', rweight := 10^-sigscores$postp[mutsig]]
    tmpo@gatk[muttype == mt & filter.reasons == 'lysis.test', rescue.fdr := 
        lysis.pv / (lysis.pv + lysis.beta * rweight * nt/na)]

    # rescue refers uniquely to rescued sites, even though regularly PASSed sites
    # would also meet these criteria.
    tmpo@gatk[muttype == mt & filter.reasons == 'lysis.test', rescue := 
        !pass & rescue.fdr <= rescue.target.fdr]
    # avoid NAs in rescue. if we really care to know that a site was also not
    # considered for rescue, we can test rescue.fdr or rweight for NA.
    tmpo@gatk[is.na(rescue), rescue := FALSE]
    data.table::setkey(tmpo@gatk, chr, pos, refnt, altnt)  # probably should already be this way

    # Now join the results back to the main (much larger) object.
    # This modifies object by reference, no need to return it.
    object@gatk[tmpo@gatk, on=.(chr, pos, refnt, altnt),
        c('rweight', 'rescue.fdr', 'rescue') := list(i.rweight, i.rescue.fdr, i.rescue)]

    # Summary info to store in the SCAN2 object's @mutsig.rescue slot.
    list(rescue.target.fdr=rescue.target.fdr,
        postp = sigscores$postp,
        test.spectrum=as.spectrum(mutsigs),
        true.sig=true.sig,
        artifact.sig=artifact.sig,
        nmuts=length(mutsigs),
        weight.true=sigscores$weight.true,
        weight.artifact=sigscores$weight.artifact,
        relative.error=sigscores$rel.error)
})



# Produces a weight vector of the same length as the provided signatures
# (which must all be of the same length). High values (>0) indicate high
# likelihood of originating from the artifact signature; low values (<0)
# indicate high likelihood of originating from the true signature.
#
# mutsigs - mutation signature channel values for each mutation under
#    consideration for rescue.  DO NOT call as.spectrum() before calling
#    get.sig.score().
get.sig.score <- function(mutsigs, true.sig, artifact.sig, eps=0.001) {
    # pracma::lsqnonneg needs a vector; a 'table' doesn't cut it
    test.spectrum <- as.numeric(as.spectrum(mutsigs))

    sigs <- cbind(true.sig, artifact.sig)
    weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
    rel.error <- norm(test.spectrum - as.vector(sigs %*% weights), type='2') / norm(test.spectrum, type='2')

    # Force weight > 0 so ratios always exist; weights can sometimes
    # be small so adding eps may have a large effect; but this is
    # somewhat controlled by using as.spectrum() on test.spectrum,
    # which converts the spectrum into a probability distn (i.e.,
    # sum(all channels) = 1).
    weights <- weights + eps
    # if we assume there are only 2 possible generating signatures,
    # then dividing weights by sum(weights) would convert the fit
    # weights into probabilities. but there's no numerical need for
    # that since we only ever consider the ratio of the weights.
    # a better future method might consider the remaining error in
    # the fit to encapsulate all other unknown processes.
    postp <- log10(artifact.sig*weights[2]) - log10(true.sig*weights[1])
    list(postp=postp, weight.true=weights[1], weight.artifact=weights[2], rel.error=rel.error)
}
