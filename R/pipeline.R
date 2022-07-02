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


# Recommended to use smaller tiles than the usual 10 MB. The files processed
# here are basepair resolution and cover essentially the entire genome.
#
# the values for the 2 following parameters are decently well tuned for the
# default n.chunks(=100) and n.permutations(=10,000).  if the ratio of
# n.permutations/n.chunks is decreased, then so should these tuning
# parameters.
# snv.N - number of random SNVs to make before each downsampling. Lower values
#         will spend more time waiting on the overhead of calls to bedtools shuffle.
# indel.K - reduces the number of random indels generated per iteration.
make.permuted.mutations <- function(sc.sample, muts, callable.bed, genome.string, genome.file, muttype=c('snv', 'indel'),
    n.permutations=10000, snv.N=1e5, indel.K=1/50, n.chunks=100, quiet=TRUE, report.mem=TRUE)
{
    muttype <- match.arg(muttype)

    cat('Permuting mutations using', n.chunks, 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    # just bins the numbers 1..n.permutations into 'n.chunks' bins, where
    # the bin size is close to equal. i.e., divvy up the number of permutations
    # to solve roughly equally across chunks.
    desired.perms <- unname(table(cut(1:n.permutations, breaks=n.chunks)))

    # Simple, not great, method for generating a sample-unique value for
    # seed construction.
    seed.base <- strtoi(paste0('0x', substr(digest::sha1(sc.sample), 1, 7)))

    progressr::with_progress({
        p <- progressr::progressor(along=1:n.chunks)
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:n.chunks, function(i) {
            pc <- perfcheck(paste('make.perms',i), {
                    if (muttype == 'snv') {
                        perms <- make.perms(muts=muts, callable=callable.bed,
                            genome.string=genome.string, genome.file=genome.file,
                            seed.base=seed.base,
                            muttype=muttype, desired.perms=desired.perms[i],
                            quiet=quiet, n.sample=snv.N)
                    } else if (muttype == 'indel') {
                        perms <- make.perms(muts=muts, callable=callable.bed,
                            genome.string=genome.string, genome.file=genome.file,
                            seed.base=seed.base,
                            muttype=muttype, desired.perms=desired.perms[i],
                            quiet=quiet, k=indel.K)
                    }
                },
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            perms
        }, future.seed=0)  # this is REQUIRED for make.perms()
    }, enable=TRUE)

    perms <- concat.perms(xs)
}


# Multicore mutation signature rescue.
#
# object.paths - character vector of paths to SCAN2 object files (.rda). IMPORTANT:
#     object.paths must be a NAMED VECTOR for which the element names point to the
#     desired OUTPUT .RDA FILES and the elements themselves are the inputs.
#
#     N.B. using compressed objects can be counterproductive. During decompression, peak
#     memory usage is size(compressed table) + size(decompressed table) + size(object).
#     Just loading a decompressed object to begin with removes the compress table size at
#     peak.  Typical values for humans are: compressed table ~ 900 MB, decompressed table ~
#     2500 MB, rest of object negligible.
# add.muts - a data.table of additional somatic mutations for creating the true somatic
#     mutation signature.  This allows the possibility of including muts from other PTA single
#     cell projects, or even different technologies (e.g., NanoSeq, META-CS), so long as the
#     added mutations are expected to share the same mutational process as the mutations
#     analyzed here.  When combining other SCAN2 runs: only "pass" mutations (i.e., VAF-based)
#     should be used; NOT SCAN2 signature-rescued mutations.
#
#     Use add.muts with care - even if the same mutational processes are active in other
#     experiments, differing technological biases or artifact processes may skew the
#     signatures.
#
#     the mutation table must contain "muttype" and "mutsig" columns.
#     Ignored if NULL.
# rescue.target.fdr - similar to the main pipeline's target.fdr. The cutoff used to rescue
#     mutations after their {lysis,mda}.fdr values have been adjusted due by mutation
#     signature rescue.
# artifact.sigs - names of signatures derived from 52
#     human neurons in Luquette et al. 2022.  Users can supply other artifact signatures
#     by providing them in the proper format (as.spectrum({sbs96,id83}(x))) in a list
#     with 'snv' and 'indel' entries.  The list elements must be the _names_ of variables
#     containing the signatures such that they can be accessed by get().
# true.sig - same format as artifact.sigs, but for the spectrum of true mutations.  Should
#     normally not be specified by the user--these are calculated from the high confidence
#     mutation calls in the objects and add.muts.
mutsig.rescue <- function(object.paths, add.muts, rescue.target.fdr=0.01,
    artifact.sigs=list(snv=data(snv.artifact.signature.v3),
                       indel=data(indel.artifact.signature.v1)),
    true.sig=NULL, quiet=FALSE, report.mem=TRUE)
{
    # Ensure that the user did set names for outputs
    if (is.null(names(object.paths)))
        stop('output RDAs must be specified in the `names()` of `object.paths`')

    # Ensure none of the output RDAs exist
    already.exists <- sapply(names(object.paths), file.exists)
    if (any(already.exists))
        stop(paste('these files specified in `names()` of `object.paths` already exist, please delete them and rerun this pipeline:', paste(names(object.paths)[already.exists], collapse=' ')))

    # Sanity check the add.muts table before doing any work
    use.add.muts <- FALSE
    if (!missing(add.muts) & !is.null(add.muts)) {
        if (!('data.table' %in% class(add.muts)) |
            !('muttype' %in% colnames(add.muts)) |
            !('mutsig' %in% colnames(add.muts)))
            stop('add.muts must be a data.table with "muttype" and "mutsig" columns')
        use.add.muts <- TRUE
    }

    cat('Rescuing mutations by signature.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')
    cat('WARNING: this pipeline requires ~6 GB of RAM per core due to data.table behavior that does not allow assignment by reference without first copying the entire table.\n')

    cat('Step 1. Getting high confidence mutations from', length(object.paths), 'SCAN2 objects.\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(object.paths))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        gatks <- future.apply::future_lapply(1:length(object.paths), function(i) {
            pc <- perfcheck(paste('prepare.object',i), {
                x <- get(load(object.paths[i]))
                if (is.compressed(x))
                    x <- decompress(x)
                ret <- reduce.table(x@gatk, target.fdr=x@call.mutations$target.fdr)
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            # Do not return anything large (like the 1-3 GB complete SCAN2 object).
            # The main R process must keep memory use low so that future.apply()
            # forked children do not accidentally copy a process with very high
            # RAM usage.
            #
            # The @gatk table is small here because prepare.object takes a small
            # subset (~1000 entries out of ~1-10 million).
            ret
        })
    }, enable=TRUE)


    cat('Step 2. Building true mutation spectra.\n')
    muttypes <- c('snv', 'indel')
    true.sigs <- setNames(lapply(muttypes, function(mt) {
        # unless user specifies it, just the raw spectrum of calls
        if (!is.null(true.sig)) {
            return(true.sig[[mt]])
        } else {
            mutsigs <- do.call(c, lapply(gatks, function(gatk)
                get.high.quality.mutations(gatk, muttype=mt)$mutsig))

            if (use.add.muts) {
                extra <- add.muts[muttype == mt]$mutsig
                cat(mt, ':', length(extra), 'mutations taken from outside sources for true signature creation (add.muts)\n')
                mutsigs <- c(mutsigs, extra)
            }

            if (mt == 'snv') true.sig <- as.spectrum(sbs96(mutsigs))
            if (mt == 'indel') true.sig <- as.spectrum(id83(mutsigs))

            cat(mt, ': created true signature from', length(mutsigs), 'high confidence mutations.\n')
            return(true.sig)
        }
    }), muttypes)

    # Slim down memory usage as much as possible before fork()ing.
    rm(gatks)
    gc()


    cat('Step3. Rescuing mutations and writing out new SCAN2 object RDA files.\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(object.paths))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        results <- future.apply::future_lapply(1:length(object.paths), function(i) {
            pc <- perfcheck(paste('mutsig.rescue.one',i), {
                x <- get(load(object.paths[i]))
                if (is.compressed(x))
                    x <- decompress(x)
                # some old objects don't have this slot; making it doesn't change the correctness of the code
                x@mutsig.rescue <- NULL   
                # there is some very odd data.table behavior that causes data.tables
                # loaded from disk to not be editable by reference.  data.table
                # recommends running setDT(), but this doesn't work for object@gatk.
                # attr(results@gatk, '.internal.selfref') is nil and it seems no
                # amount of messing with it will allow me to edit it.
                # I have also tried setalloccol() with no luck; changing mutsig.rescue.one() to
                # return a small data.table which is then joined onto x@gatk in this function
                # to avoid issues with function scope; etc.  Nothing worked except the below
                # (which is VERY bad): using <- assignment (an equivalent := assignment by
                # reference does NOT work).
                x@gatk$rescue <- FALSE   # this will be updated by mutsig.rescue.one where appropriate
                # Why is this bad? it copies the data.table, doubling the memory use
                # of a 2.5-3.0 GB object. This simply won't be usable on something like
                # a mouse SCAN2 object with full SNP table.
                #
                # The unfortunate result of all this is we need ~6 GB per core to
                # comfortably run this pipeline for human cells.
                #
                # use options(datatable.verbose = TRUE) for debugging
    
                for (mt in muttypes) {
                    x@mutsig.rescue[[mt]] <- mutsig.rescue.one(x,
                        muttype=mt,
                        artifact.sig=get(artifact.sigs[[mt]]),
                        true.sig=true.sigs[[mt]],
                        rescue.target.fdr=rescue.target.fdr)
                }
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            results <- x
            save(results, file=names(object.paths)[i], compress=FALSE)

            # add in sample ID and bulk ID in anticipation of the batch-wide
            # mutation tables.
            results@gatk$sample <- results@single.cell
            results@gatk$bulk.sample <- results@bulk
            # can't have sample names in columns - rbind would like uniform column names
            colnames(results@gatk)[colnames(results@gatk) == results@single.cell] <- 'scgt'
            ret <- list(
                sample=results@single.cell,
                muts=results@gatk[pass == TRUE | rescue == TRUE],
                # one entry for each mutation type (snv, indel)
                sig.homogeneity.test=
                    lapply(results@mutsig.rescue, function(msr) msr$sig.homogeneity.test)
            )
            ret
        })
    }, enable=TRUE)

    list(muts=do.call(rbind, lapply(results, function(r) r$muts)),
         sig.homogeneity.tests=setNames(lapply(results, function(r) r$sig.homogeneity.test),
                                        sapply(results, function(r) r$sample)))
}



# NOT concat.perms - this function combines permutations across samples, not
# across parallelized chunks of permutation-generation.
combine.permutations <- function(perm.files, genome.string, report.mem=TRUE) {
    genome.seqinfo <- genome.string.to.seqinfo.object(genome.string)

    pc <- perfcheck('load permutations', {
        data <- future.apply::future_lapply(perm.files, function(f) {
            load(f)  # loads permdata
            #cat(length(permdata$perms), 'permutations\n')
            list(
                seeds.used=data.frame(sample=permdata$muts[1,]$sample,
                           file=f, seed.used=permdata$seeds.used),
                # keep only necessary columns to reduce mem usage
                perms=lapply(permdata$perms, function(p) p[, c('chr', 'pos', 'mutsig')])
            )
        })
    }, report.mem=report.mem)
    cat(pc, '\n')
    
    seeds.used <- do.call(rbind, lapply(data, function(d) d$seeds.used))
    if (sum(duplicated(seeds.used)) != 0)
        warning('detected duplicate seeds, there may be a problem in how you supplied seed.base to permtool.R!')

    perml <- lapply(data, function(d) d$perms)
    if (any(sapply(perml, length) != length(perml[[1]])))
        stop('all permutation files must contain the same number of permutations')

    # perml can be very large (~5G for the 52 PTA paper neurons) and needs to
    # be copied once for every thread. There is probably a better way to divide
    # this to avoid copies.
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(perml[[1]]))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        zperml <- GenomicRanges::GRangesList(future.apply::future_lapply(1:length(perml[[1]]), function(i) {
            pc <- perfcheck(paste('restructure permutations', i), {
                z <- lapply(perml, function(ps) ps[[i]])
                names(z) <- NULL # GRanges c() won't combine things with different names
                # the permutations are now dataframes, not GRanges, so rbind and convert
                z <- do.call(rbind, z)
                gz <- GenomicRanges::GRanges(seqnames=z$chr, ranges=IRanges(start=z$pos, width=1),
                    seqinfo=genome.seqinfo)
                gz <- sort(GenomeInfoDb::sortSeqlevels(gz))
                gz$mutsig <- z$mutsig
                gz$perm.id <- i
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            gz
        }))
    })
    list(seeds.used=seeds.used, zperml=zperml)
}
