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
    setNames(lapply(object.paths, function(path) {
        if (!quiet) print(path)
        x <- get(load(path))
        prepare.object(x)
    }), sapply(os, function(o) o@single.cell))
}


# object.paths - character vector of paths to SCAN2 object files (.rda)
# objects - list of SCAN2 objects produced by get.objects.for.sig.rescue(). This
#     can take a long time to load, so users may wish to call get.objects.for.sig.rescue()
#     on their own.
# rescue.target.fdr - similar to the main pipeline's target.fdr. The cutoff used to rescue
#     mutations after their {lysis,mda}.fdr values have been adjusted due by mutation
#     signature rescue.
# artifact.sigs - names of signatures derived from 52
#     human neurons in Luquette et al. 2022.  Users can supply other artifact signatures
#     by providing them in the proper format (as.spectrum({sbs96,id83}(x))) in a list
#     with 'snv' and 'indel' entries.  The list elements must be the _names_ of variables
#     containing the signatures such that they can be accessed by get().
# true.sig - the spectrum of VAF-based mutation calls in the SCAN2 objects in object.paths.
#     Can also be overridden if desired.
mutsig.rescue <- function(object.paths, rescue.target.fdr=0.01, objects=NULL,
    artifact.sigs=list(snv=data(snv.artifact.signature.v3), indel=data(indel.artifact.signature.v1)),
    true.sig=NULL, quiet=FALSE)
{
    if (!missing(object.paths) & !is.null(objects))
        stop('exactly one of object.paths or objects may be specified')

    # Both of these code paths create copies of the SCAN2 objects with a
    # very small subset of the @gatk table selected.
    if (is.null(objects)) {
        objects <- get.objects.for.sig.rescue(object.paths, quiet=quiet)
    } else {
        objects <- setNames(lapply(objects, prepare.object, quiet=quiet),
            sapply(objects, function(o) o@single.cell))
    }

    if (!quiet)
        cat(sprintf('SCAN2 signature-based rescue on %d samples\n', length(objects)))

    muttypes <- c('snv', 'indel')
    summaries <- setNames(lapply(muttypes, function(mt) {
        artifact.sig <- get(artifact.sigs[[mt]])

        # unless user specifies it, just the raw spectrum of calls
        if (is.null(true.sig)) {
            mutsigs <- do.call(c, lapply(objects, function(o) o@gatk[muttype == mt & pass == TRUE]$mutsig))
            if (mt == 'snv') true.sig <- as.spectrum(sbs96(mutsigs))
            if (mt == 'indel') true.sig <- as.spectrum(id83(mutsigs))

            cat(mt, ': created true signature from', length(mutsigs), 'high confidence mutations.\n')
        }

        cat(sprintf('%s : adjusting FDR for lysis artifacts using new FDR target=%0.4f..\n',
            mt, rescue.target.fdr))

        ret <- setNames(lapply(objects, function(o) {
            cat('', o@single.cell)
            mutsig.rescue.one(o, muttype=mt,
                artifact.sig=artifact.sig, true.sig=true.sig, rescue.target.fdr=rescue.target.fdr)
        }), sapply(objects, function(o) o@single.cell))
        cat('.\n')
        ret
    }), muttypes)

    # All work done above modified objects by reference, now recollating summary data
    # This will be added as a SCAN2 slot one day.
    lapply(objects, function(o) { list(object=o, snv=summaries[['snv']][[o@single.cell]], indel=summaries[['indel']][[o@single.cell]]) })
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




# XXX: in the future it'd be nice to record mutation signature rescue and
# summaries in a new SCAN2 object slot, but the issue is RAM.  SCAN2
# objects (for humans) are large (2.5-3 GB) and we need to read through
# a good number of them (i.e., 30-100 for various projects) to construct
# the true mutation signature.  Storing these all in memory for later isn't
# possible and re-reading them all can take an hour, so for now mutation
# rescue ifno just isn't stored in the object.
setGeneric("mutsig.rescue.one", function(object, artifact.sig, true.sig, rescue.target.fdr=0.01, muttype=c('snv', 'indel'), ...)
    standardGeneric("mutsig.rescue.one"))
setMethod("mutsig.rescue.one", "SCAN2", function(object, artifact.sig, true.sig, rescue.target.fdr=0.01, muttype=c('snv', 'indel'), ...) {
    mt <- match.arg(muttype)

    sigtype <- if (mt == 'snv') sbs96 else id83
    mutsigs <- sigtype(object@gatk[muttype == mt & filter.reasons == 'lysis.test']$mutsig)

    sigscores <- get.sig.score(mutsigs=mutsigs,
        artifact.sig=artifact.sig, true.sig=true.sig, ...)

    # it doesn't seem to be possible to use a column assigned by := for another
    # assignment in the same data.table statement.  i.e., to combine all of these
    # into a single statement.
    object@gatk[muttype == mt & filter.reasons == 'lysis.test', rweight := 10^-sigscores$postp[mutsig]]
    object@gatk[muttype == mt & filter.reasons == 'lysis.test', rescue.fdr := 
        lysis.pv / (lysis.pv + lysis.beta * rweight * nt/na)]

    # rescue refers uniquely to rescued sites, even though regularly PASSed sites
    # would also meet these criteria.
    object@gatk[muttype == mt & filter.reasons == 'lysis.test', rescue := 
        !pass & rescue.fdr <= rescue.target.fdr]

    # no need to return object since all changes were made by reference
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
