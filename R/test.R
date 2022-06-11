# some utility functions for comparing results
check.length <- function(a, b, msg) {
    ret <- length(a) != length(b)
    if (ret)
        cat(sprintf('FAILED lengths (a=%d, b=%d): %s', length(a), length(b), msg))
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
        w <- which(xor(is.na(a), is.na(b)) | abs(a - b) > 
            tolerance)
        cat('indexes: '); print(w)
        cat('abs(diffs): '); print(abs(a-b)[w])
    } else cat('.')
}


testpipe <- function(test.data=c('legacy_tiny', 'legacy_chr22', 'legacy_custom'), verbose=FALSE, custom=NULL, legacy=TRUE, n.cores=future::availableCores())
{
    if (n.cores > 1) {
        require(future)
        future::plan(multicore, workers=n.cores)
    }

    test.data <- match.arg(test.data)
    if (test.data == 'legacy_custom' & is.null(custom))
        stop('when test.data=legacy_custom, a list must be supplied to "custom"')

    # test parameters and files
    if (test.data != 'legacy_custom') {
        sc.sample <- 'h25'
        bulk.sample <- 'hunamp'
        fpath <- function(...) system.file('extdata', paste0(test.data, '_', ...), package='scan2')
    } else {
        sc.sample <- custom$sc.sample
        bulk.sample <- custom$bulk.sample
        fpath <- function(...) paste0(custom$path, '/', ...)
    }

    mmq60 <- fpath('mmq60.tab.gz')
    mmq1 <- fpath('mmq1.tab.gz')
    phased.vcf <- fpath('phased_all.vcf.gz')
    panel <- fpath('cross_sample_panel.tab.gz')
    abfits <- fpath('fits.rda')
    sccigars <- fpath('sc_somatic_and_hsnp_spikein_cigars.tab.gz')
    bulkcigars <- fpath('bulk_somatic_and_hsnp_spikein_cigars.tab.gz')
    trainingcigars <- fpath('cigardata.tab.gz')

    if (test.data == 'legacy_tiny') {
        grs <- GRanges(seqnames=22, ranges=IRanges(start=c(30e6, 31e6),
            end=c(30999999, 31999999)))
    } else if (test.data == 'legacy_chr22') {
        # Cover chr22 with 1mb tiles
        grs <- tileGenome(seqlengths=seqinfo(genome.string.to.bsgenome.object('hs37d5'))[as.character(22)], tilewidth=1e6, cut.last.tile.in.chrom=TRUE)
    } else if (test.data == 'legacy_custom') {
        grs <- custom$grs
    }

    int.tab.path <- tempfile()
    int.tab.gz.path <- paste0(int.tab.path, '.gz')
    x <- make.integrated.table(mmq60.tab=mmq60, mmq1.tab=mmq1, phased.vcf=phased.vcf,
        bulk.sample=bulk.sample, genome='hs37d5', grs=grs)
    write.integrated.table(inttab=x$gatk, out.tab=int.tab.path, out.tab.gz=int.tab.gz.path)

    run.pipeline(sc.sample=sc.sample, bulk.sample=bulk.sample,
        int.tab=int.tab.gz.path, abfits=abfits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        trainingcigars=trainingcigars, panel=panel,
        legacy=legacy, genome='hs37d5', grs=grs, verbose=FALSE)
}

test.output <- function(pipeline.output, test.data=c('legacy_tiny', 'legacy_chr22', 'legacy_custom'), custom=NULL) {
    test.data <- match.arg(test.data)

    for (mt in c('snv', 'indel')) {
        if (test.data != 'legacy_custom') {
            legacy.rda <-
                system.file('extdata', paste0(test.data, '_', mt, '_somatic_genotypes.rda'),
                    package='scan2')
        } else {
            legacy.rda <- paste0(custom$path, '/', mt, '_somatic_genotypes.rda')
        }

        l <- get(load(legacy.rda))
        cat(' MUTTYPE =', mt, '-------------------------------------------\n')

        # necessary to look at a subset of the legacy output and pipeline
        # output. in the new pipeline's legacy mode, sites that are not admitted
        # for FDR prior estimation are *also* not scored for FDR. these sites are
        # are filtered out by both pipelines in the final calling due to static
        # filters so ignoring them will not affect validity. the sites can't be
        # compared due to missing FDR scores in the new pipeline.
        #
        # don't need to use the higher DP >= 10 cutoff for indels because that
        # was applied after processing, so the lower DP sites will be present.
        l <- l[l$bulk.dp >= 11 & (is.na(l$alt.1.lowmq) | l$alt.1.lowmq == 0) & l$dp >= 6,]
    
        # always use snv filters here: legacy code did not apply DP >= 10 (the
        # only difference between SNV and indel hard filters) until after analysis.
        sfp <- pipeline.output@static.filter.params[['snv']]
        p <- pipeline.output@gatk[
                muttype == mt &
                (muttype == 'snv' | !is.na(nalleles)) & # sites not in the panel were removed by merge(), but only for indels
                bulk.dp >= sfp$min.bulk.dp &
                (is.na(balt.lowmq) | balt.lowmq == 0) &
                balt == 0 & bulk.gt == '0/0' &
                (dbsnp == '.' | !sfp$exclude.dbsnp) &
                scalt >= sfp$min.sc.alt &
                dp >= sfp$min.sc.dp]

        test.equal(l$chr, p$chr, "chr")
        test.equal(l$pos, p$pos, "pos")
        test.equal(l$refnt, p$refnt, "refnt")
        test.equal(l$altnt, p$altnt, "altnt")
        test.equal(l$dbsnp, p$dbsnp, "dbsnp")
        test.equal(l$h25, p$h25, "h25")
        test.equal(l$hunamp, p$bulk.gt, "bulk.gt")  # used to be called hunamp, now unambiguously labeled as bulk.gt
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

        # test cross sample panel values
        if (mt == 'indel') {
            test.equal(l$nalleles, p$nalleles, 'panel: nalleles')
            test.equal(l$unique.donors, p$unique.donors, 'panel: unique.donors')
            test.equal(l$unique.cells, p$unique.cells, 'panel: unique.cells')
            test.equal(l$unique.bulks, p$unique.bulks, 'panel: unique.bulks')
            test.equal(l$max.out, p$max.out, 'panel: max.out')
            test.equal(l$sum.out, p$sum.out, 'panel: sum.out')
            test.equal(l$sum.bulk, p$sum.bulk, 'panel: sum.bulk')
        }
        # CIGARs are necessarily different because the legacy script (which used
        # samtools view at every candidate site and was too slow for all-sites mode)
        # produces different CIGAR counts than the new script (which uses pysam).
        # the counts generally trend together very well, but they would have to be
        # exact for these tests to work out.
        if (FALSE) {
            test.tol(l$id.score.y, p$id.score.y, "id.score.y")
            test.tol(l$id.score.x, p$id.score.x, "id.score.x")
            test.tol(l$id.score, p$id.score, "id.score")
            test.tol(l$hs.score.y, p$hs.score.y, "hs.score.y")
            test.tol(l$hs.score.x, p$hs.score.x, "hs.score.x")
            test.tol(l$hs.score, p$hs.score, "hs.score")
            test.tol(l$cigar.id.test, p$cigar.id.test, "cigar.id.test")
            test.tol(l$cigar.hs.test, p$cigar.hs.test, "cigar.hs.test")
        
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
        }
        test.equal(l$lowmq.test, p$lowmq.test, "lowmq.test")
        test.equal(l$dp.test, p$dp.test, "dp.test")
        cat('\n')
    }
}
