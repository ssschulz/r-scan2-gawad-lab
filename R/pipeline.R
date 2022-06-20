# Run `expr` and then print a few statistics about memory and runtime
# print.header - if TRUE, msg and expr are ignored.  It does NOT mean
#     print msg/expr and ADD a header!
# report.mem - don't run gc(). each gc() takes 1-2s; on normal whole
#     genome workloads this is unnoticable, but when running small test
#     pipelines this can make up 99%+ of the runtime.
perfcheck <- function(msg, expr, print.header=FALSE, report.mem=TRUE) {
    if (print.header) {
        return(sprintf('%30s | %9s %11s %9s %9s',
            'Step (chunk)', 'Mem Mb', 'Peak mem Mb', 'Time (s)', 'Elapsed'))
    }
    t <- system.time(eval(expr), gcFirst=report.mem)
    mem.used <- NA
    max.mem.used <- NA
    if (report.mem) {
        g <- gc(full=TRUE, reset=TRUE)
        mem.used <- sum(g[,which(colnames(g)=='used')+1])
        max.mem.used <- sum(g[,which(colnames(g)=='max used')+1])
    }
    sprintf('%30s |  %7.1f %10.1f %7.1f %7.1f', msg,
        mem.used, max.mem.used,
        # combine user, system, and child cpu time
        sum(t[names(t) != 'elapsed']),
        t['elapsed'])
}


run.pipeline <- function(
    sc.sample, bulk.sample,
    int.tab,
    abfits,
    sccigars, bulkcigars, trainingcigars,
    dptab,
    genome,
    genome.seqinfo=genome.string.to.seqinfo.object(genome),
    config.yaml=NULL,
    grs=tileGenome(seqlengths=genome.seqinfo[as.character(1:22)], tilewidth=10e6, cut.last.tile.in.chrom=TRUE),
    legacy=FALSE, report.mem=TRUE, verbose=TRUE)
{
    cat('Starting chunked SCAN2 pipeline on', length(grs), 'chunks\n')
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores\n')
    cat('Detailed chunk schedule:\n')
    cat(sprintf('%7s %5s %10s %10s\n', 'Chunk', 'Chr', 'Start', 'End'))
    for (i in 1:length(grs)) {
        cat(sprintf('%7d %5s %10d %10d\n',
            i, as.character(seqnames(grs)[i]), start(grs)[i], end(grs)[i]))
    }

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]
            # Don't put the perfcheck() calls in p(), because progressr
            # doesn't evaluate those arguments if progress bars re turned off.
            pc <- perfcheck(paste('make.scan',i),
                gt <- make.scan(single.cell=sc.sample, bulk=bulk.sample, genome=genome, region=gr), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            if (!is.null(config.yaml)) {
                gt <- add.static.filter.params(gt, config.path=config.yaml)
            } else {
                # the default values in this function ARE NOT GUARANTEED
                # to match the defaults in the scan2 script.
                warning('config.yaml not specified, using R function default static filter parameters, which may not match the scan2 tool')
                gt <- add.static.filter.params(gt, muttype='snv')
                gt <- add.static.filter.params(gt, muttype='indel')
            }

            pc <- perfcheck(paste('read.integrated.table',i),
                gt <- read.integrated.table(gt, path=int.tab, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.ab.fits',i),
                gt <- add.ab.fits(gt, path=abfits), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.ab.estimates',i),
                gt <- compute.ab.estimates(gt, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.cigar.data',i),
                gt <- add.cigar.data(gt, sccigars, bulkcigars, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.models',i),
                gt <- compute.models(gt, verbose=verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            # Note about legacy mode: legacy mode isn't actually legacy, it's what
            # is still currently used. Legacy uses only resampled germline hets for
            # the null CIGAR distn while the new mode uses all germline hets. The new
            # mode is *way* too slow to ever use because the computation is O(n) where
            # n is the number of null sites. The new mode needs to approximate the
            # 2-d CIGAR op probability space with a fixed N (e.g., of gaussians) to
            # guarantee reasonable runtime.
            pc <- perfcheck(paste('compute.excess.cigar.scores',i),
                gt <- compute.excess.cigar.scores(object=gt, path=trainingcigars, legacy=TRUE, quiet=!verbose),
                    report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.static.filters',i),
                gt <- compute.static.filters(gt), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            p()
            gt
        }, future.seed=0)  # CRITICAL! library(future) ensures that each child process
                           # has a different random seed.
    })
    cat("Chunked pipeline complete.\n")

    x <- do.call(concat, xs)
    x <- call.mutations(x, target.fdr=0.01, mode=ifelse(legacy, 'legacy', 'new'), quiet=!verbose)
    x <- add.depth.profile(x, depth.path=dptab)
    x <- compute.mutburden(x)
    x
}



# Full GATK annotation pipeline. Creates an annotated integrated table, which
# contains many site-specific annotations and the full matrix of alt and ref
# read counts for all single cells and bulks.
make.integrated.table <- function(mmq60.tab, mmq1.tab, phased.vcf,
    bulk.sample, genome, genome.seqinfo=genome.string.to.seqinfo.object(genome), panel=NULL,
    grs=tileGenome(seqlengths=genome.seqinfo[as.character(1:22)], tilewidth=10e6, cut.last.tile.in.chrom=TRUE),
    quiet=TRUE, report.mem=FALSE)
{
    cat('Starting integrated table pipeline on', length(grs), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]

            pc <- perfcheck(paste('read and annotate raw data',i), {
                gatk <- read.tabix.data(path=mmq60.tab, region=gr, quiet=quiet,
                    colClasses=list(character='chr'))  # force chromosome column to be type=str

                # Columns in GATK are split as site data | sample-specific count data/genotypes
                # There are 7 site data columns (chr, pos, dbsnp ID, ref allele, alt allele, mq, mqrs).
                # Try to keep the columns split by site-wide data | sample-specific data
                sitewide <- gatk[,1:7]
                samplespecific <- gatk[,-(1:7)]

                annotate.gatk.bulk(sitewide, samplespecific, bulk.sample, quiet=quiet)
                annotate.gatk(gatk=sitewide, gatk.counts=samplespecific, genome.string=genome, add.mutsig=TRUE)
                annotate.gatk.lowmq(sitewide, path=mmq1.tab, bulk=bulk.sample, region=gr, quiet=quiet)
                annotate.gatk.phasing(sitewide, phasing.path=phased.vcf, region=gr, quiet=quiet)
                annotate.gatk.panel(sitewide, panel.path=panel, region=gr, quiet=quiet)
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            cbind(sitewide, samplespecific)
        })
    })

    gatk <- rbindlist(xs)
    resampling.details <- gatk.resample.phased.sites(gatk)
    list(gatk=gatk, resampling.details=resampling.details)
}


digest.depth.profile <- function(path, sc.sample, bulk.sample,
    genome, genome.seqinfo=genome.string.to.seqinfo.object(genome),
    clamp.dp=500,
    grs=tileGenome(seqlengths=genome.seqinfo[as.character(1:22)], tilewidth=10e6, cut.last.tile.in.chrom=TRUE),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Digesting depth profile using', length(grs), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]

            pc <- perfcheck(paste('digest.depth.2sample',i),
                dptab <- digest.depth.2sample(path=path, sc.sample=sc.sample,
                    bulk.sample=bulk.sample, clamp.dp=clamp.dp, region=gr, quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            dptab
        })
    }, enable=TRUE)

    # Sum all of the tables
    dptab <- Reduce(`+`, xs)

    list(dptab=dptab, clamp.dp=clamp.dp)
}


# Recommended to use smaller tiles than the usual 10 MB. The files processed
# here are basepair resolution and cover essentially the entire genome.
make.callable.regions <- function(path, sc.sample, bulk.sample,
    genome, min.sc.dp, min.bulk.dp,
    genome.seqinfo=genome.string.to.seqinfo.object(genome),
    grs=tileGenome(seqlengths=genome.seqinfo[as.character(1:22)], tilewidth=5e6, cut.last.tile.in.chrom=TRUE),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Getting callable regions using', length(grs), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]

            pc <- perfcheck(paste('compute.callable.region',i),
                g <- compute.callable.region(path=path, sc.sample=sc.sample,
                    bulk.sample=bulk.sample, min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp,
                    region=gr, quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            g
        })
    }, enable=TRUE)

    gr <- do.call(c, xs)
    list(regions=gr, sc.sample=sc.sample, bulk.sample=bulk.sample, min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp)
}
