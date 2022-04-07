testpipe <- function(test.data=c('legacy_tiny', 'legacy_chr22', 'legacy_custom'), verbose=FALSE, custom=NULL) {
    test.data <- match.arg(test.data)
    if (test.data == 'legacy_custom' & is.null(custom))
        stop('when test.data=legacy_custom, a data.frame must be supplied to "custom"')

    # test parameters and files
    if (legacy.data != 'legacy_custom') {
        sc.sample <- 'h25'
        bulk.sample <- 'hunamp'
        fpath <- function(...) system.file('extdata', paste0(test.data, '_', ...), package='scan2')
    } else {
        sc.sample <- custom$sc.sample
        bulk.sample <- custom$bulk.sample
        fpath <- function(...) paste0(custom$path, '/', ...)
    }

    mmq60 <- fpath('mmq60.tab.bgz"')
    mmq1 <- fpath('mmq1.tab.bgz')
    hsnps <- fpath('hsnps.tab.bgz')
    abfits <- fpath('fits.rda')
    sccigars <- fpath('h25_somatic_and_hsnp_spikein_cigars.tab.bgz')
    bulkcigars <- fpath('hunamp_somatic_and_hsnp_spikein_cigars.tab.bgz')


    require(future)
    require(future.apply)
    plan(multisession)
    #plan(sequential) # just for testing
    require(progressr)
    if (test.data == 'legacy_tiny') {
        grs <- list(GRanges(seqnames=22, ranges=IRanges(start=30e6, end=30999999)),
                    GRanges(seqnames=22, ranges=IRanges(start=31e6, end=31999999)))
    } else if (test.data == 'legacy_chr22') {
        grs <- lapply(16:49, function(start)
            GRanges(seqnames=22, ranges=IRanges(start=start*1e6, (start+1)*1e6-1)))
    } else if (test.data == 'legacy_custom') {
        grs <- custom$grs
    }

    printfun <- invisible
    if (verbose)
        printfun <- print

    perfcheck <- function(msg, expr) {
        t <- system.time(eval(expr))
        g <- gc(reset=TRUE)
        sprintf('%30s |  %10.1f %10.1f %9.2f', msg,
            sum(g[,which(colnames(g)=='used')+1]),
            sum(g[,which(colnames(g)=='max used')+1]),
            # combine user, system, and child cpu time
            sum(t[names(t) != 'elapsed']))
    }

    cat('Starting chunked pipeline on', length(grs), 'chunks\n')
    cat('Parallelizing with', future::availableCores(), 'cores\n')
    cat('Detailed chunk schedule:\n')
    cat(sprintf('%7s %5s %10s %10s\n', 'Chunk', 'Chr', 'Start', 'End'))
    for (i in 1:length(grs)) {
        cat(sprintf('%7d %5s %10d %10d\n', i, seqnames(grs[[i]])[1],
            start(grs[[i]]), end(grs[[i]])))
    }

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky',
            sprintf('%30s |  %11s %11s %11s', 'Step (chunk)', 'Mem Mb', 'Peak mem Mb', 'Time s'))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[[i]]
            p(class='sticky', amount=0, perfcheck(paste('make.scan',i),
                gt <- make.scan(single.cell=sc.sample, bulk=bulk.sample, genome='hs37d5', region=gr)))
            gt <- add.static.filter.params(gt)
            p(class='sticky', amount=0, perfcheck(paste('read.gatk',i),
                x1 <- read.gatk(gt, path=mmq60, quiet=!verbose)))
            p(class='sticky', amount=0, perfcheck(paste('read.gatk.lowmq',i),
                y1 <- read.gatk.lowmq(x1, path=mmq1, quiet=!verbose)))
            p(class='sticky', amount=0, perfcheck(paste('add.training.data',i),
                z1 <- add.training.data(y1, path=hsnps, quiet=!verbose)))
            p(class='sticky', amount=0, perfcheck(paste('add.ab.fits',i),
                w1 <- add.ab.fits(z1, path=abfits)))
            p(class='sticky', amount=0, perfcheck(paste('compute.ab.estimates',i),
                v1 <- compute.ab.estimates(w1, quiet=!verbose)))
            p(class='sticky', amount=0, perfcheck(paste('add.cigar.data',i),
                r1 <- add.cigar.data(v1, sccigars, bulkcigars, quiet=!verbose)))
            p(class='sticky', amount=0, perfcheck(paste('compute.models',i),
                s1 <- compute.models(r1)))
            p()
            s1
        })
    })

    x <- do.call(concat, xs)
    printfun(system.time(x2 <- resample.training.data(x)))
    printfun('After resample.training.data:\n') ; printfun(gc(reset=TRUE))
    printfun(system.time(x4 <- compute.excess.cigar.scores(x2)))
    printfun('After compute.excess.cigar.scores:\n') ; printfun(gc(reset=TRUE))
    printfun(system.time(x5 <- compute.static.filters(x4)))
    printfun('After compute.static.filters:\n') ; printfun(gc(reset=TRUE))
    printfun(system.time(x6 <- compute.fdr.priors(x5)))
    printfun('After compute.fdr.priors:\n') ; printfun(gc(reset=TRUE))
    printfun(system.time(x7 <- compute.fdr(x6)))
    printfun('After compute.fdr:\n') ; printfun(gc(reset=TRUE))

    list(xs, x7)
}


test.output <- function(pipeline.output, custom.path, test.data=c('legacy_tiny', 'legacy_chr22', 'legacy_custom')) {
    test.data <- match.arg(test.data)

    if (test.data != 'legacy_custom') {
        legacy.rda <-
            system.file('extdata', paste0(test.data, '_somatic_genotypes.rda'),
                package='scan2')
    } else {
        legacy.rda <- paste0(custom.path, '/', 'somatic_genotypes.rda')
    }
    l <- get(load(legacy.rda, verb=F))

    # necessary to look at a subset of the legacy output and pipeline
    # output. in the new pipeline's legacy mode, sites that are not admitted
    # for FDR prior estimation are *also* not scored for FDR. these sites are
    # are filtered out by both pipelines in the final calling due to static
    # filters so ignoring them will not affect validity. the sites can't be
    # compared due to missing FDR scores in the new pipeline.
    l <- l[l$bulk.dp >= 11 & (is.na(l$alt.1.lowmq) | l$alt.1.lowmq == 0) & l$dp >= 6,]
    p <- pipeline.output@gatk[!is.na(lysis.fdr)]

    check.length <- function(a, b, msg) {
        ret <- length(a) != length(b)
        if (ret)
            cat(paste('FAILED lengths (a=%d, b=%d): %s', length(a), length(b), msg))
        return(ret)
    }

    test.equal <- function(a, b, msg) {
        if (check.length(a, b, msg))
            return()
        nfail <- sum(xor(is.na(a), is.na(b)) | a != b, na.rm=TRUE)
        if (nfail > 0)
            cat(paste0('FAILED (',nfail,') equality:', msg, '\n'))
        else cat('.')
    }
    test.tol <- function(a, b, msg, tolerance=1e-6) {
        if (check.length(a, b, msg))
            return()
        nfail <- sum(xor(is.na(a), is.na(b)) | abs(a-b) > tolerance, na.rm=TRUE)
        if (nfail > 0) {
            cat(sprintf('FAILED (%d) tolerance (%f):', nfail, tolerance), msg, '\n')
            print(which(
                xor(is.na(a), is.na(b)) | abs(a-b) > tolerance))
        } else cat('.')
    }

    test.equal(l$chr, p$chr, "chr")
    test.equal(l$pos, p$pos, "pos")
    test.equal(l$refnt, p$refnt, "refnt")
    test.equal(l$altnt, p$altnt, "altnt")
    test.equal(l$dbsnp, p$dbsnp, "dbsnp")
    test.equal(l$h25, p$h25, "h25")
    test.equal(l$hunamp, p$hunamp, "hunamp")
    test.equal(l$dp, p$dp, "dp")
    test.equal(l$af, p$af, "af")
    test.equal(l$bulk.dp, p$bulk.dp, "bulk.dp")

    test.tol(abs(l$gp.mu), abs(p$gp.mu), "gp.mu")
    test.tol(l$gp.sd, p$gp.sd, "gp.sd")
    test.tol(pmin(l$ab,1-l$ab), pmin(p$ab,1-p$ab), "ab")
    test.tol(l$abc.pv, p$abc.pv, "abc.pv")
    test.tol(l$lysis.pv, p$lysis.pv, "lysis.pv")
    test.tol(l$mda.pv, p$mda.pv, "mda.pv")
    test.tol(l$nt, p$nt, "nt")
    test.tol(l$na, p$na, "na")
    # the beta in the current SCAN2 table is not derived from the same min.
    # FDR method used in legacy. the old beta is computed internally when
    # calculating the final FDR estimates, but it is not saved in the table.
    #test.tol(l$lysis.beta, p$lysis.beta, "lysis.beta")
    test.tol(l$lysis.fdr, p$lysis.fdr, "lysis.fdr")
    #test.tol(l$mda.beta, p$mda.beta, "mda.beta")
    test.tol(l$mda.fdr, p$mda.fdr, "mda.fdr")
    test.tol(l$id.score.y, p$id.score.y, "id.score.y")
    test.tol(l$id.score.x, p$id.score.x, "id.score.x")
    test.tol(l$id.score, p$id.score, "id.score")
    test.tol(l$hs.score.y, p$hs.score.y, "hs.score.y")
    test.tol(l$hs.score.x, p$hs.score.x, "hs.score.x")
    test.tol(l$hs.score, p$hs.score, "hs.score")
    test.tol(l$cigar.id.test, p$cigar.id.test, "cigar.id.test")
    test.tol(l$cigar.hs.test, p$cigar.hs.test, "cigar.hs.test")

    test.equal(l$lowmq.test, p$lowmq.test, "lowmq.test")
    test.equal(l$M.cigars, p$M.cigars, "M.cigars")
    test.equal(l$ID.cigars, p$ID.cigars, "ID.cigars")
    test.equal(l$HS.cigars, p$HS.cigars, "HS.cigars")
    test.equal(l$other.cigars, p$other.cigars, "other.cigars")
    test.equal(l$dp.cigars, p$dp.cigars, "dp.cigars")
    test.equal(l$M.cigars.bulk, p$M.cigars.bulk, "M.cigars.bulk")
    test.equal(l$ID.cigars.bulk, p$ID.cigars.bulk, "ID.cigars.bulk")
    test.equal(l$HS.cigars.bulk, p$HS.cigars.bulk, "HS.cigars.bulk")
    test.equal(l$other.cigars.bulk, p$other.cigars.bulk, "other.cigars.bulk")
    test.equal(l$dp.cigars.bulk, p$dp.cigars.bulk, "dp.cigars.bulk")
    test.equal(l$dp.test, p$dp.test, "dp.test")
    cat('\n')
}
