run.pipeline <- function(
    sc.sample, bulk.sample,
    mmq60, mmq1,
    hsnps,
    abfits,
    sccigars, bulkcigars,
    genome, grs=tileGenome(seqlengths=seqinfo(genome.string.to.bsgenome.object(genome))[as.character(1:22)], tilewidth=10e6, cut.last.tile.in.chrom=TRUE),
    verbose=TRUE)
{
    cat('Starting chunked SCAN2 pipeline on', length(grs), 'chunks\n')
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)
    cat('Parallelizing with', future::availableCores(), 'cores\n')

    printfun <- invisible
    if (verbose)
        printfun <- print

    perfcheck <- function(msg, expr) {
        t <- system.time(eval(expr))
        g <- gc(reset=TRUE)
        sprintf('%30s |  %7.1f %10.1f %7.1f %7.1f', msg,
            sum(g[,which(colnames(g)=='used')+1]),
            sum(g[,which(colnames(g)=='max used')+1]),
            # combine user, system, and child cpu time
            sum(t[names(t) != 'elapsed']),
            t['elapsed'])
    }

    cat('Detailed chunk schedule:\n')
    cat(sprintf('%7s %5s %10s %10s\n', 'Chunk', 'Chr', 'Start', 'End'))
    for (i in 1:length(grs)) {
        cat(sprintf('%7d %5s %10d %10d\n',
            i, as.character(seqnames(grs)[i]), start(grs)[i], end(grs)[i]))
    }

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky',
            sprintf('%30s | %9s %11s %9s %9s', 'Step (chunk)', 'Mem Mb', 'Peak mem Mb', 'Time (s)', 'Elapsed'))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]
            pc <- perfcheck(paste('make.scan',i),
                gt <- make.scan(single.cell=sc.sample, bulk=bulk.sample, genome=genome, region=gr))
            p(class='sticky', amount=0, pc)

            gt <- add.static.filter.params(gt)

            pc <- perfcheck(paste('read.gatk',i),
                x1 <- read.gatk(gt, path=mmq60, quiet=TRUE))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('read.gatk.lowmq',i),
                y1 <- read.gatk.lowmq(x1, path=mmq1, quiet=TRUE))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.training.data',i),
                z1 <- add.training.data(y1, path=hsnps, quiet=!verbose))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.ab.fits',i),
                w1 <- add.ab.fits(z1, path=abfits))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.ab.estimates',i),
                v1 <- compute.ab.estimates(w1, quiet=!verbose))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.cigar.data',i),
                r1 <- add.cigar.data(v1, sccigars, bulkcigars, quiet=!verbose))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.models',i),
                s1 <- compute.models(r1))
            p(class='sticky', amount=0, pc)
            p()
            s1
        }, future.seed=0)  # CRITICAL! library(future) ensures that each child process
                           # has a different random seed.
    })

    x <- do.call(concat, xs)
    perfcheck('resample.training.data', x2 <- resample.training.data(x))
    perfcheck('compute.excess.cigar.scores', x4 <- compute.excess.cigar.scores(x2))
    perfcheck('compute.static.filters', x5 <- compute.static.filters(x4))
    perfcheck('compute.fdr.priors', x6 <- compute.fdr.priors(x5))
    perfcheck('compute.fdr', x7 <- compute.fdr(x6))

    list(xs, x7)
}
