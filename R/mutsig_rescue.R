# Steps:
#   1. add filter reasons to each sample (only add reasons for sites
#      with >=min.alt alt reads because it's slow.)
#   2. get the "good signature"
#        this should really be handled by NMF
#   3. run rescue()
do.rescue <- function(ldf, lysis.sig, good.sig, min.alt=2, calling.fdr=0.01, rescue.fdr=0.01) {
    # step 1
    cat(sprintf('SCAN2 multi-sample rescue on %d samples, %d high confidence SNVs from target.fdr=%0.4f\n',
        length(ldf), sum(sapply(ldf, function(d) sum(d$pass))), calling.fdr))

    cat('    Step 1: annotating filter reasons and SBS96 status..\n       ')
    ldf2 <- lapply(1:length(ldf), function(i) {
        cat(sprintf(' %s', names(ldf)[i]))
        df <- ldf[[i]]
        # reduce the size of the dataframe before adding reasons, major speed up
        df <- df[!is.na(df$af) & round(df$af*df$dp) >= min.alt,]
        df$filter.reasons <- get.filter.reasons(df,
            min.alt=min.alt,
            fdr.threshold=calling.fdr)
        get.3mer(df)
    })
    names(ldf2) <- names(ldf)
    ldf <- ldf2
    cat('.\n')

    # step 2
    # unless user specifies it, just the raw spectrum of calls
    cat('    Step 2: constructing true SNV spectrum..\n')
    if (missing(good.sig)) {
        x <- do.call(rbind, lapply(ldf, function(s) s[s$pass,1:5]))
        good.sig <- df.to.sbs96(x)
    }

    # step 3
    cat(sprintf('    Step 3: adjusting FDR for lysis artifacts using new FDR target=%0.4f..\n        ',
        rescue.fdr))
    final <- lapply(1:length(ldf), function(i) {
        df <- ldf[[i]]
        cat(sprintf(' %s', names(ldf)[i]))
        rescue(df, lysis.sig=lysis.sig, good.sig=good.sig, rescue.fdr=rescue.fdr)
    })
    names(final) <- names(ldf)
    cat('.\n')

    list(calling.fdr=calling.fdr,                    # for posterity
        rescue.fdr=rescue.fdr,                       #
        good.sig = good.sig, lysis.sig = lysis.sig,  #
        df = lapply(final, function(l) l$df),
        postp = lapply(final, function(l) l$postp),
        nsnvs = sapply(final, function(l) l$nsnvs),
        weight.true = sapply(final, function(l) l$weight.true),
        weight.artifact = sapply(final, function(l) l$weight.artifact))
}



# alt.idx is usually either 10 or 13, depending on whether the bulk
# sample name is lexicographically before the single cell sample name
# must set:
#   - fdr.threshold = calling threshold
#   - min.alt = calling threshold
get.filter.reasons <- function(df, fdr.threshold=0.01, min.alt=2) {
    # the alt column isn't in a consistent place, so recover it from VAF
    alt.reads <- round(df$af*df$dp)
    m <- cbind(
        pass=df$pass,
        abc.test=df$abc.pv <= 0.05,
        lysis.test=df$lysis.fdr > fdr.threshold,
        mda.test=df$mda.fdr > fdr.threshold,
        cigar.ID=!df$cigar.id.test,
        cigar.HS=!df$cigar.hs.test,
        lowmq=!df$lowmq.test,
        dp=!df$dp.test,
        min.alt=is.na(alt.reads) | alt.reads < min.alt # NA implies DP=0 implies nalt=0
    )
    apply(m, 1, function(row) paste(colnames(m)[row], collapse='&'))
}



df.to.sbs96 <- function(df) {
    if (!('type.and.ctx' %in% colnames(df)))
        df <- get.3mer(df)

    # copied from plot.3mer
    bases <- c("A", "C", "G", "T")
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each = 4), rep(c("C", "T"), 
        each = 48), rep(bases, times = 4), ":", rep(c("C", "T"), 
        each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
        rep(c("A", "C", "G"), each = 16)))
    t2 <- table(df$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[, 2])]
    t <- t + 0.1 # add a pseudocount for 0s 
    t/sum(t)
}



rescue <- function(df, lysis.sig, good.sig, rescue.fdr, ...) {
    sigscores <- get.sig.score(ssnvs=df[df$filter.reason=='lysis.test',],
        lysis.sig=lysis.sig, good.sig=good.sig, ...)
    postp <- sigscores$postp
    df$rweight <- 10^-postp[df$type.and.ctx]
    df$lysis.fdr2 <- df$lysis.alpha /
        (df$lysis.alpha + df$lysis.beta * df$rweight * df$nt/df$na)
    # pass2 refers uniquely to rescued sites
    df$pass2 <- !df$pass & df$lysis.fdr2 <= rescue.fdr & df$filter.reasons == 'lysis.test'
    df$filter.reasons[df$pass2] <-
        paste0(df$filter.reasons[df$pass2], ';rescue')
    list(df = df, postp = postp,
        nsnvs=sigscores$nsnvs,
        weight.true=sigscores$weight.true,
        weight.artifact=sigscores$weight.artifact)
}



get.sig.score <- function(ssnvs, good.sig, lysis.sig, eps=0.001) {
    test.sig <- df.to.sbs96(ssnvs)
    sigs <- cbind(good.sig, lysis.sig)
    weights <- pracma::lsqnonneg(sigs, test.sig)$x
    recon <- as.vector(sigs %*% weights)

    # Make sure no channels = 0
    weights <- weights + eps
    weights <- weights / sum(weights)
    nsnvs <- sum(ssnvs$filter.reason=='lysis.test')

    postp <- log10(lysis.sig*weights[2]) - log10(good.sig*weights[1])
    list(postp=postp, nsnvs=nsnvs, weight.true=weights[1], weight.artifact=weights[2])
}
