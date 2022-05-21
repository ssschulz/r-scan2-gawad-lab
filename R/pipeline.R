# Run `expr` and then print a few statistics about memory and runtime
# print.header - if TRUE, msg and expr are ignored.  It does NOT mean
#     print msg/expr and ADD a header!
perfcheck <- function(msg, expr, print.header=FALSE) {
    if (print.header) {
        return(sprintf('%30s | %9s %11s %9s %9s',
            'Step (chunk)', 'Mem Mb', 'Peak mem Mb', 'Time (s)', 'Elapsed'))
    }
    t <- system.time(eval(expr))
    g <- gc(reset=TRUE)
    sprintf('%30s |  %7.1f %10.1f %7.1f %7.1f', msg,
        sum(g[,which(colnames(g)=='used')+1]),
        sum(g[,which(colnames(g)=='max used')+1]),
        # combine user, system, and child cpu time
        sum(t[names(t) != 'elapsed']),
        t['elapsed'])
}


run.pipeline <- function(
    sc.sample, bulk.sample,
    mmq60, mmq1,
    hsnps,
    abfits,
    sccigars, bulkcigars, trainingcigars,
    fdr.prior.data,
    genome,
    tmpsave.rda,
    grs=tileGenome(seqlengths=seqinfo(genome.string.to.bsgenome.object(genome))[as.character(1:22)], tilewidth=10e6, cut.last.tile.in.chrom=TRUE),
    verbose=TRUE)
{
    if (!missing(tmpsave.rda) & file.exists(tmpsave.rda))
        stop('temporary save file tmpsave.rda already exists, please delete it first')

    cat('Starting chunked SCAN2 pipeline on', length(grs), 'chunks\n')
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)
    cat('Parallelizing with', future::availableCores(), 'cores\n')
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
                z1 <- add.training.data(y1, path=hsnps, quiet=!verbose, require.resampled=TRUE))
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
                s1 <- compute.models(r1, verbose=verbose))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.excess.cigar.scores',i),
                t1 <- compute.excess.cigar.scores(s1, trainingcigars, quiet=!verbose))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.static.filters',i),
                u1 <- compute.static.filters(t1))
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.fdr',i),
                v1 <- compute.fdr(u1, fdr.prior.data, mode='new'))
            p(class='sticky', amount=0, pc)

            p()
            v1
        }, future.seed=0)  # CRITICAL! library(future) ensures that each child process
                           # has a different random seed.
    })
    cat("Chunked pipeline complete.\n")
    if (!missing(tmpsave.rda)) {
        cat('Saving pre-merged chunks to', tmpsave.rda, '\n')
        save(xs, file=tmpsave.rda)
    }

    x <- do.call(concat, xs)
    cat("Merged SCAN2 object after chunked pipeline:\n")
    print(x)

    x
}
