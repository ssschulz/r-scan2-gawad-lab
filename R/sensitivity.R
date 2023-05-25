# tilewidth purposefully not exposed here
compute.hsnp.spatial.sens <- function(object, muttype=c('snv', 'indel'), quiet=TRUE) {
    muttype <- match.arg(muttype)

    if (!quiet) cat("Tiling genome with 10 kb tiles\n")
    tiles <- tileGenome(seqlengths=seqlengths(results@region), tilewidth=1e3, cut.last.tile.in.chrom=T)

    # For hSNP sensitivity with 10kb windows, acf() suggests useful info only
    # exists for the nearest tile up- or downstream. Even so, the correlation
    # is only 0.1.
    compute.hsnp.sens.for.tiles(object, tiles, smooth.tiles=1, muttype=muttype, quiet=quiet)
}


compute.hsnp.sens.for.tiles <- function(object, tiles, smooth.tiles=0, muttype=c('snv', 'indel'), quiet=TRUE) {
    muttype <- match.arg(muttype)

    # neither muttype == muttype nor muttype == ..muttype work. whatever.
    different.variable.name <- muttype
    hsnps <- object@gatk[training.site == TRUE & muttype == different.variable.name]
    gr <- GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
        seqinfo=object@genome.seqinfo)
    gr$pass <- hsnps$training.pass

    if (!quiet) cat(paste0("Mapping ", nrow(hsnps), " germline ", muttype, "s to ", length(tiles), " non-overlapping tiles\n"))
    tiles$n.called <- countOverlaps(tiles, gr[gr$pass == TRUE])
    tiles$n.tested <- countOverlaps(tiles, gr)
    # this will be 0, just to allow for smoothing, then set to NA
    tiles$sens <- ifelse(tiles$n.tested > 0, tiles$n.called/tiles$n.tested, 0)
    tiles$sens.smoothed <- as.numeric(stats::filter(tiles$sens,
        filter=c(1, rep(1, 2*smooth.tiles))/(1+2*smooth.tiles), sides=2))
    tiles$sens <- ifelse(tiles$n.tested > 0, tiles$sens, NA)
    tiles
}



# Spatial sensitivity applies only to VAF-based calling! Not to mutation signature-based rescue!
compute.spatial.sensitivity.depth <- function(single.cell.id, bulk.id,
    static.filter.params, joint.dptab.path, genome.string,
    grs.for.sens=genome.string.to.tiling(genome.string, tilewidth=1e3, group='auto'),
    grs.for.parallelization=genome.string.to.tiling(genome.string, tilewidth=10e6, group='auto'),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Gathering read depth data for spatial somatic calling sensitivity using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    sfp <- static.filter.params

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            # Only retain grs.for.sens windows that *start* in the
            # grs.for.parallelization window. This avoids assigning a grs.for.sens
            # window to 2 grs.for.parallelization windows if it spans the boundary.
            # There's probably a better way to do this.
            grs.for.sens2 <- GRanges(seqnames=seqnames(grs.for.sens),
                ranges=IRanges(start=start(grs.for.sens), width=1),
                seqinfo=seqinfo(grs.for.sens))
            gr2 <- grs.for.sens[IRanges::countOverlaps(grs.for.sens2, gr, minoverlap=1) > 0,]

            pc <- perfcheck(paste('read.depth.2sample', i),
                dp <- read.depth.2sample(path=joint.dptab.path, sc.sample=single.cell.id,
                    bulk.sample=bulk.id, keep.coords=TRUE, region=reduce(gr2), quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', pc, amount=0)

            # Set up several GRanges objects for the binnedAverage() method
            pc <- perfcheck(paste("compute.averages", i), {
                    dp.gr <- GRanges(seqnames=dp$chr, ranges=IRanges(dp$pos, width=1),
                        sc.dp=dp[[3]],
                        bulk.dp=dp[[4]],
                        # these are logical vectors: 1 if each single base is >=
                        # the filter cutoff. The binned average can then be multiplied by the width
                        # to recover number of bases.
                        base.gt.snv.sc.min.dp=dp[[3]] >= sfp$snv$min.sc.dp,
                        base.gt.snv.bulk.min.dp=dp[[4]] >= sfp$snv$min.bulk.dp,
                        base.gt.indel.sc.min.dp=dp[[3]] >= sfp$indel$min.sc.dp,
                        base.gt.indel.bulk.min.dp=dp[[4]] >= sfp$indel$min.bulk.dp,
                        seqinfo=seqinfo(gr2))
                    sc.dp <- mcolAsRleList(dp.gr, 'sc.dp')
                    bulk.dp <- mcolAsRleList(dp.gr, 'bulk.dp')
                    base.gt.snv.sc.min.dp <- mcolAsRleList(dp.gr, 'base.gt.snv.sc.min.dp')
                    base.gt.snv.bulk.min.dp <- mcolAsRleList(dp.gr, 'base.gt.snv.bulk.min.dp')
                    base.gt.indel.sc.min.dp <- mcolAsRleList(dp.gr, 'base.gt.indel.sc.min.dp')
                    base.gt.indel.bulk.min.dp <- mcolAsRleList(dp.gr, 'base.gt.indel.bulk.min.dp')
                    gr2 <- binnedAverage(gr2, sc.dp, 'mean.sc.dp')
                    gr2 <- binnedAverage(gr2, bulk.dp, 'mean.bulk.dp')
                    gr2 <- binnedAverage(gr2, base.gt.snv.sc.min.dp, 'mean.gt.snv.sc.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.snv.bulk.min.dp, 'mean.gt.snv.bulk.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.indel.sc.min.dp, 'mean.gt.indel.sc.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.indel.bulk.min.dp, 'mean.gt.indel.bulk.min.dp')
                }, report.mem=report.mem)
            p(class='sticky', pc, amount=1)

            data.table(chr=as.character(seqnames(gr2)), start=start(gr2), end=end(gr2),
                mean.sc.dp=gr2$mean.sc.dp, mean.bulk.dp=gr2$mean.bulk.dp,
                bases.gt.snv.sc.min.dp=width(gr2) * gr2$mean.gt.snv.sc.min.dp,
                bases.gt.snv.bulk.min.dp=width(gr2) * gr2$mean.gt.snv.bulk.min.dp,
                bases.gt.indel.sc.min.dp=width(gr2) * gr2$mean.gt.indel.sc.min.dp,
                bases.gt.indel.bulk.min.dp=width(gr2) * gr2$mean.gt.indel.bulk.min.dp
            )
        })
    }, enable=TRUE)

    return(rbindlist(xs))  # decide if we want to keep this as a GRanges object or not
}



# Spatial sensitivity applies only to VAF-based calling! Not to mutation signature-based rescue!
# IMPORTANT: grs.for.sens must match the GRanges used for compute.spatial.sensitivity.depth().
# N.B. do NOT try to rewrite this function to use a SCAN2 object as argument. future()'s multi-
# core implementation copies the entire parent memory space to each child thread no matter how I
# try to prevent it. For human genomes, this means a waste of ~2-2.5G of RAM per thread.
compute.spatial.sensitivity.abmodel <- function(
    single.cell.id, ab.fits, integrated.table.path, genome.string,
    grs.for.sens=genome.string.to.tiling(genome.string, tilewidth=1e3, group='auto'),
    grs.for.parallelization=genome.string.to.tiling(genome.string, tilewidth=10e6, group='auto'),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Estimating genome-wide AB for spatial somatic calling sensitivity using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')
    # Necessary for AB estimation (infer.gp uses matrix multiplication)
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)

    cat('Estimating AB at window mid-points..\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            # Only retain grs.for.sens windows that *start* in the
            # grs.for.parallelization window. This avoids assigning a grs.for.sens
            # window to 2 grs.for.parallelization windows if it spans the boundary.
            # There's probably a better way to do this.
            grs.for.sens2 <- GRanges(seqnames=seqnames(grs.for.sens),
                ranges=IRanges(start=start(grs.for.sens), width=1),
                seqinfo=seqinfo(grs.for.sens))
            gr2 <- grs.for.sens[IRanges::countOverlaps(grs.for.sens2, gr, minoverlap=1) > 0,]

            pc <- perfcheck(paste('get.training.sites',i),
                    training.sites <- get.training.sites.for.abmodel.by.range(
                        region=reduce(gr2), integrated.table.path=integrated.table.path,
                        single.cell.id=single.cell.id, quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.ab',i),
                    ab <- compute.ab.given.sites.and.training.data(
                        sites=data.table(chr=as.character(seqnames(gr2)), pos=(end(gr2)+start(gr2)) / 2),
                        training.hsnps=training.sites,
                        ab.fits=ab.fits, quiet=TRUE),  # have to set quiet=TRUE or progress bar will get overridden
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            cbind(data.table(chr=as.character(seqnames(gr2)), start=start(gr2), end=end(gr2)), ab)
        })
    }, enable=TRUE)

    return(rbindlist(xs))
}



# abmodel.covs - compute.spatial.sensitivity.abmodel() output
# depth.covs - compute.spatial.sensitivity.depth() output
# This function doesn't use future() for multicore support, so passing the large
# SCAN2 `object` is feasible.
integrate.spatial.sensitivity.covariates <- function(object, abmodel.covs, depth.covs,
    grs.for.sens=genome.string.to.tiling(object@genome.string, tilewidth=1e3, group='auto'))
{
    cat('Counting germline training hSNP sites..\n')
    ab.data <- data.table(chr=as.character(seqnames(grs.for.sens)),
        start=start(grs.for.sens), end=end(grs.for.sens))
    hsnps <- object@gatk[training.site == TRUE & muttype == 'snv']
    ab.data$n.training.hsnps <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    hsnps <- object@gatk[training.site == TRUE & muttype == 'snv' & af >= 0.5]
    ab.data$n.training.hsnps.maj <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    hsnps <- object@gatk[training.site == TRUE & muttype == 'snv' & af < 0.5]
    ab.data$n.training.hsnps.min <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    # now just the called sites, for sens estimation
    hsnps <- object@gatk[training.pass == TRUE & muttype == 'snv']
    ab.data$n.training.hsnps.passed <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    hsnps <- object@gatk[training.pass == TRUE & muttype == 'snv' & af >= 0.5]
    ab.data$n.training.hsnps.maj.passed <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    hsnps <- object@gatk[training.pass == TRUE & muttype == 'snv' & af < 0.5]
    ab.data$n.training.hsnps.min.passed <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hsnps$chr, ranges=IRanges(start=hsnps$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))

    cat('Counting germline het indel sites..\n')
    hindels <- object@gatk[training.site == TRUE & muttype == 'indel']
    ab.data$n.training.hindels <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hindels$chr, ranges=IRanges(start=hindels$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))
    hindels <- hindels[training.pass == TRUE]  # now just the called sites, for sens estimation
    ab.data$n.training.hindels.passed <- IRanges::countOverlaps(grs.for.sens,
        GRanges(seqnames=hindels$chr, ranges=IRanges(start=hindels$pos, width=1),
            seqinfo=genome.string.to.seqinfo.object(object@genome.string)))

    cat("Merging with AB model covariates..\n")
    ret <- merge(ab.data, abmodel.covs, by=c('chr', 'start', 'end'))
    cat("Merging with depth covariates..\n")
    ret <- merge(ret, depth.covs, by=c('chr', 'start', 'end'))
    ret
}
