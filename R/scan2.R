setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClassUnion('null.or.BSgenome', c('NULL', 'BSgenome'))
setClassUnion('null.or.GRanges', c('NULL', 'GRanges'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.logical', c('NULL', 'logical'))
setClassUnion('null.or.character', c('NULL', 'character'))

# if adding new slots to the class, be sure to update concat(...)
# to combine the slots properly.
setClass("SCAN2", slots=c(
    region='null.or.GRanges',
    genome.string='character',
    genome.object='null.or.BSgenome',
    single.cell='character',
    bulk='character',
    gatk="null.or.df",
    integrated.table.path='null.or.character',
    resampled.training.data="null.or.list",
    ab.fits='null.or.df',
    static.filter.params='null.or.list',
    ab.estimates='null.or.df',
    mut.models='null.or.df',
    cigar.data='null.or.df',
    excess.cigar.scores='null.or.df',
    fdr.prior.data='null.or.list',
    fdr='null.or.list'))


# "flip" the GP mu around 0 to best match the AF of each candidate mutation. 
match.ab <- function(af, gp.mu) {
    ifelse(af < 1/2, -abs(gp.mu), abs(gp.mu))
}

# currently will always fire when SCAN2 is built from chunked objects.
# firing isn't a problem, just causes a warning message.
# need to define what the "full genome" is (in terms of GRanges objects)
# for the recognized genomes.
check.chunked <- function(object, message) {
    if (!is.null(object@region))
        warning(message)
}

genome.string.to.bsgenome.object <- function(genome=c('hs37d5', 'hg38', 'mm10')) {
    if (genome == 'hs37d5') {
        require(BSgenome.Hsapiens.1000genomes.hs37d5)
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else if (genome == 'hg38') {
        require(BSgenome.Hsapiens.NCBI.GRCh38)
        genome <- BSgenome.Hsapiens.NCBI.GRCh38
    } else if (genome == 'mm9') {
        require(BSgenome.Mmusculus.UCSC.mm10)
        genome <- BSgenome.Mmusculus.UCSC.mm10
    } else {
        stop(paste('genome string', genome, 'is not supported'))
    }
    genome
}


make.scan <- function(single.cell, bulk, genome=c('hs37d5', 'hg38', 'mm10'), region=NULL) {
    genome <- match.arg(genome)
    new("SCAN2", single.cell=single.cell, bulk=bulk,
        genome.string=genome,
        genome.object=genome.string.to.bsgenome.object(genome),
        region=region,
        gatk=NULL,
        ab.fits=NULL,
        # these slots used to hold tables; now their data is incorporated into @gatk and
        # they only record analysis parameters, if any apply.
        ab.estimates=NULL,
        mut.models=NULL,
        cigar.data=NULL,
        excess.cigar.scores=NULL,
        static.filter.params=NULL,
        fdr.prior.data=NULL,
        fdr=NULL)
}


setValidity("SCAN2", function(object) {
    if (length(object@genome.string) != 1)
        return("must provide exactly one genome name")

    if (!is.null(object@genome.object)) {
        if (class(object@genome.object) != 'BSgenome')
            stop('@genome.object must be of class BSgenome')
    }

    if (length(object@single.cell) != 1)
        return("must provide exactly one single cell sample name")

    if (length(object@bulk) != 1)
        return("must provide exactly one bulk sample name")

    if (!is.null(object@region)) {
        if (!is(object@region, 'GRanges')) 
            stop('"region" must be a GRanges object')
        # This isn't possible if we're to combine chunks for the final
        # pipeline steps. Need a check that allows either single region
        # GRanges or the expected full autosome set.
        #
        # Disabling for now, but this isn't safe in general.
        #if (length(object@region) != 1) {
           #stop('"region" must contain exactly 1 range')
        #}
    }

    if (!is.null(object@gatk)) {
        if (ncol(object@gatk) != 13) {
            return("@gatk must have 13 columns")
        }
        if (!all(colnames(object@gatk)[1:7] ==
                c('chr', 'pos', 'refnt', 'altnt', 'dbsnp', 'mq', 'mqrs'))) {
            return("@gatk is improperly formatted")
        }
    }

    if (!is.null(object@ab.fits)) {
        # XXX: should also check that relevant chromosomes are present
        if (!all(colnames(object@ab.fits == c('a', 'b', 'c', 'd', 'logp'))))
            return("@ab.fits is improperly formatted")
    }

    if (!is.null(object@ab.estimates)) {
        if (!all(c('ab', 'gp.mu', 'gp.sd') %in% colnames(object@gatk)))
            return("ab.estimates were not properly imported")
    }

    if (!is.null(object@mut.models)) {
        if (!all(c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')
            %in% colnames(object@gatk)))
            return("mut.models were not properly calculated")
    }

    # XXX: add validator for rest of slots
    return(TRUE)
})


# various steps in the pipeline require certain slots to be filled.
# this function will ensure all slots in 'slots: (character)' are not
# NULL and will generate an error message otherwise.
check.slots <- function(object, slots, abort=TRUE) {
    error.occurred = FALSE
    for (s in slots) {
        if (is.null(slot(object, s))) {
            error.occurred = TRUE
            if (s == 'gatk')
                cat("must import GATK read counts first (see: read.gatk())\n")
            if (s == 'ab.fits')
                cat('must import AB model parameters (see: add.ab.fits())\n')
            if (s == 'ab.estimates')
                cat("must import allele balance model estimates first (see: add.ab.estimates())\n")
            if (s == 'cigar.data')
                cat("must import CIGAR data first (see: add.cigar.data())\n")
            if (s == 'excess.cigar.scores')
                cat("must compute excess CIGAR scores first (see: compute.excess.cigar.scores())\n")
            if (s == 'static.filter' | s == 'static.filters' | s == 'static.filter.params')
                cat("must apply static site filters first (see: add.static.filters())\n")
            if (s == 'fdr.prior.data')
                cat("must compute or import FDR priors first (see: compute.fdr.priors())\n")
        }
    }
    if (error.occurred & abort)
        stop("One or more required slots are missing. See above for details.")
    return(error.occurred)
}


setMethod("show", "SCAN2", function(object) {
    cat("#", is(object)[[1]], "\n")
    if (!is.null(object@region)) {
        cat("#   Region:")
        if (length(object@region) > 1) {
            cat("\n")
            print(object@region)
        } else {
            cat('',length(object@region),'intervals\n')
        }
    }
    cat("#   Genome:", object@genome.string, "\n")
    cat("#   Single cell ID:", object@single.cell, "\n")
    cat("#   Bulk ID:", object@bulk, "\n")
    cat("#   GATK:")
    if (is.null(object@gatk)) {
        cat(" (no data)\n")
    } else {
        cat('', nrow(object@gatk), "raw sites\n")
    }

    cat("#   AB model training hSNPs:")
    if (!('training.site' %in% object@gatk)) {
        cat(" (no data)\n")
    } else {
        # germline indels are not used for AB model training
        tdata <- object@gatk[training.site==TRUE & muttype=='snv']
        cat('', nrow(tdata),
            sprintf("phased sites (hap1=%d, hap2=%d)",
                sum(tdata$phased.gt=='1|0', na.rm=TRUE),
                sum(tdata$phased.gt=='0|1', na.rm=TRUE)),"\n")
        zzz <- do.call(rbind, lapply(split(tdata, tdata$chr),
            function(td) {
                af <- td$phased.hap1/(td$phased.hap1+td$phased.hap2)
                cbind(dist=diff(td$pos), af1=af[-nrow(td)], af2=td$af[-1])
        }))
        cors <- sapply(10^(2:5), function(threshold) {
            z <- zzz[zzz[,1] <= threshold,]; cor(z[,2], z[,3], use='complete.obs') }) 
        cat('#       VAF correlation between neighboring hSNPs:\n')
        cat('#           <100 bp', round(cors[1], 3),
            '<1000 bp', round(cors[2], 3),
            '<10 kbp', round(cors[3], 3),
            '<100 kbp', round(cors[4], 3), '\n')
        if ('resampled.training.site' %in% colnames(object@gatk)) {
            cat('#        ', nrow(object@gatk[resampled.training.site == TRUE & muttype == 'snv']),
                'resampled hSNPs\n')
        }
    }

    cat("#   AB model parameters:")
    if (is.null(object@ab.fits)) {
        cat(" (no data)\n")
    } else {
        cat(sprintf('\n#       average (over chromosomes): a=%0.3f, b=%0.3f, c=%0.3f, d=%0.3f\n',
            mean(object@ab.fits$a),
            mean(object@ab.fits$b),
            mean(object@ab.fits$c),
            mean(object@ab.fits$d)))
    }

    cat("#   Allele balance:")
    if (is.null(object@ab.estimates)) {
        cat(" (not computed)\n")
    } else {
        s <- summary(object@gatk$gp.sd)
        cat('\n#       mean (0 is neutral):',
            round(mean(object@gatk$gp.mu), 3), '\n')
        cat('#       uncertainty (Q25, median, Q75):',
            round(s['1st Qu.'], 3),
            round(s['Median'], 3),
            round(s['3rd Qu.'], 3), '\n')
        if ('training.site' %in% colnames(object@gatk)) {
            xs <- round(object@gatk[training.site==TRUE & muttype == 'snv',
                .(mean=mean(gp.mu), cor=cor(af, ab, use='complete.obs'))],3)
            cat('#       mean at training hSNPs:', xs$mean, '\n')
            # computing correlation doesn't make sense here because AF and AB
            # do not necessarily refer to the same haplotype: AF is always the
            # mutated haplotype.
                #.(mean=mean(gp.mu), cor=cor(af, ab, use='complete.obs'))],3)
            #cat('#       correlation with VAF at training hSNPs', xs$cor, '\n')
        }
    }

    cat("#   Mutation models:")
    if (is.null(object@mut.models)) {
        cat(" (not computed)\n")
    } else {
        cat(" computed\n")
    }

    cat("#   CIGAR data:")
    if (is.null(object@cigar.data)) {
        cat(" (no data)\n")
    } else {
        cat(' single cell:', object@cigar.data$sc.sites, "sites, ")
        cat('bulk:', object@cigar.data$bulk.sites, "sites\n")
    }

    cat("#   Static filters: ")
    if (!('static.filter' %in% colnames(object@gatk))) {
        cat("(not applied)\n")
    } else {
        cat(sum(object@gatk$static.filter, na.rm=TRUE), "retained",
            sum(!object@gatk$static.filter, na.rm=TRUE), "removed",
            sum(is.na(object@gatk$static.filter)), "NA\n")
    }
})


# to use: e.g., x <- do.call(concat, xs)
# where xs is a list of SCAN2 chunks. the chunks must not overlap. one day this
# might be enforced programatically
setGeneric("concat", function(...) standardGeneric("concat"))
setMethod("concat", signature="SCAN2", function(...) {
    args <- list(...)

    if (length(args) == 0)
        stop('tried to concat() 0 objects')

    if (length(args) == 1)
        return(args[[1]])

    init <- args[[1]]

    ret <- make.scan(init@single.cell, init@bulk,
        genome=init@genome.string,
        region=reduce(do.call(c, lapply(args, function(a) a@region))))

    # rbindlist quickly concatenates data.tables
    ret@gatk <- rbindlist(lapply(args, function(a) a@gatk))

    ensure.same <- function(l, slot.name, var.name) {
        # either all slots are NULL or they have the same value
        if (all(sapply(l, function(element) is.null(slot(element, slot.name))))) {
            return()
        }
        # if var.name isn't supplied, assume that `slot.name` is just a vector
        if (!missing(var.name)) {
            if (!all(sapply(l, function(element)
                all(slot(l[[1]], slot.name)[[var.name]] == slot(element, slot.name)[[var.name]]))))
                stop(paste('list of SCAN2 chunks cannot be concatenated; slot', slot.name,
                    'does not have consistent values for', var.name))
        } else {
            if (!all(sapply(l, function(element)
                all(slot(l[[1]], slot.name) == slot(element, slot.name)))))
                stop(paste('list of SCAN2 chunks cannot be concatenated; slot', slot.name,
                    'does not have consistent values'))
        }
    }

    # XXX: TODO: would be nice to ensure same genome object. just not sure how to
    # check equality at the moment.
    ensure.same(args, 'genome.string')
    ensure.same(args, 'single.cell')
    ensure.same(args, 'bulk')
    ensure.same(args, 'static.filter.params', 'min.sc.alt')
    ensure.same(args, 'static.filter.params', 'min.sc.dp')
    ensure.same(args, 'static.filter.params', 'min.bulk.dp')
    ensure.same(args, 'static.filter.params', 'max.bulk.alt')
    ensure.same(args, 'static.filter.params', 'exclude.dbsnp')
    ensure.same(args, 'static.filter.params', 'cg.id.q')
    ensure.same(args, 'static.filter.params', 'cg.hs.q')
    ret@static.filter.params <- init@static.filter.params
    ensure.same(args, 'cigar.data', 'sc.path')
    ensure.same(args, 'cigar.data', 'bulk.path')
    ret@cigar.data <- data.frame(sc.sites=sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$sc.sites))),
                               bulk.sites=sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$bulk.sites)))
    )
    ensure.same(args, 'excess.cigar.scores', 'legacy')
    ensure.same(args, 'excess.cigar.scores', 'training.sites')
    ensure.same(args, 'fdr.prior.data', 'bins')
    ensure.same(args, 'fdr.prior.data', 'max.dp')
    ensure.same(args, 'fdr.prior.data', 'candidates.used')
    ensure.same(args, 'fdr.prior.data', 'hsnps.used')
    ensure.same(args, 'fdr.prior.data', 'nt.tab')
    ensure.same(args, 'fdr.prior.data', 'na.tab')
    ensure.same(args, 'fdr.prior.data', 'mode')
    # fdr.prior.data: 'fcs' should also be identical, but the list is a little
    # inconvenient to check. the above should detect misuse 99% of the time.
    ret@fdr.prior.data <- init@fdr.prior.data
    ensure.same(args, 'fdr', 'mode')
    ret@fdr <- data.frame(mode=init@fdr$mode,
        sites=sum(sapply(args, function(a) ifelse(is.null(a@fdr), 0, a@fdr$sites))))

    # policy: ab.fits has to be the same for all chunks being concat()ed
    # if that's true, just use the first one
    if (any(sapply(args, function(a) any(as.matrix(init@ab.fits) != as.matrix(a@ab.fits)))))
        stop('@ab.fits must be identical for all concat() elements')
    ret@ab.fits <- init@ab.fits

    ret@ab.estimates <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@ab.estimates), 0, a@ab.estimates$sites))))

    ret@mut.models <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@mut.models), 0, a@mut.models$sites))))

    ret
})



setGeneric("read.integrated.table", function(object, path, quiet=FALSE)
    standardGeneric("read.integrated.table"))
setMethod("read.integrated.table", "SCAN2", function(object, path, quiet=FALSE) {
    object@gatk <- read.and.annotate.integrated.table(path=path, sample.id=object@single.cell,
        region=object@region, quiet=quiet)
    object@integrated.table.path <- path
    object
})



# Add AB Gaussian process parameter fits for each chromosome.
# These take a long time to compute, so the snakemake pipeline will almost
# always be necessary here.
setGeneric("add.ab.fits", function(object, path)
    standardGeneric("add.ab.fits"))
setMethod("add.ab.fits", "SCAN2", function(object, path) {
    object@ab.fits <- get(load(path))
    object
})



# Using tile subsampling, AB model fitting can feasibly be performed on
# a single multicore node.
# AB fits must be computed using all training sites on a chromosome.
# AB fits are produced by sampling 20,000 random parameter 4-tuples. The
# highest logP values are considered the best fits. The random sampling
# is iteratively refined in 'n.steps' steps by restricting random parameter
# sampling to a smaller subset of the space.
#   - samples.per.chunk, n.chunks: compute 'samples.per.chunk' random
#         parameter samplings in 'n.chunks' independent (possibly parallel)
#         threads. samples.per.chunk*n.chunks should be kept at 20,000.
#   - refine.n.steps: number of iterations in which the (a,b,c,d) parameter space
#         is refined and sampled.
#   - refine.top.n: use the top 'top.n' parameter values (by logP) to create the
#         next refined parameter space
#   - n.tiles and hsnp.tilesize: the log-likelihood of (a,b,c,d|training hSNPs)
#         is approximated by breaking the hSNPs into non-overlapping tiles of
#         size hsnp.tilesize. This is necessary because the approximation requires
#         inverting the (#hSNPs x #hSNPs) covariance matrix. Furthermore, the
#         additional information about (a,b,c,d) provided by each tile of hSNPs
#         becomes less and less as more tiles are added. We have found that ~200
#         tiles (=20,000 hSNPs) is a good trade-off between compute time and
#         accuracy.
#   - alim, blim, clim, dlim - starting bounds for the (a,b,c,d) parameter
#         space.
setGeneric("compute.ab.fits", function(object, path, chroms=1:22,
    n.cores=future::availableCores(),
    logp.samples.per.step=20000, refine.n.steps=4, refine.top.n=50,
    n.tiles=250, hsnp.tilesize=100,
    alim=c(-7, 2), blim=c(2, 4), clim=c(-7, 2), dlim=c(2, 6))
    standardGeneric("compute.ab.fits"))
setMethod("compute.ab.fits", "SCAN2", function(object, path, chroms=1:22,
    n.cores=future::availableCores(),
    logp.samples.per.step=20000, refine.n.steps=4, refine.top.n=50,
    n.tiles=250, hsnp.tilesize=100,
    alim=c(-7, 2), blim=c(2, 4), clim=c(-7, 2), dlim=c(2, 6)) 
{
    cat("using", n.cores, "cores\n")
    n.chunks <- 100  # 4*n.cores  # using a multiple of n.cores gives a smoother progress bar
    # using n.cores and ceiling can lead to (very small) differences in number
    # of samples taken, and thus slightly different results. just use a very large
    # number of chunks (like 100)
    if (logp.samples.per.step %% n.chunks != 0)
        stop(sprintf('logp.samples.per.step must be a multiple of n.chunks (%d)', n.chunks))
    samples.per.chunk <- ceiling(logp.samples.per.step / n.chunks)

    # Just for convenience. Allow "chroms=1:22" to work for autosomes
    # Check all chroms up front so the loop doesn't die after a significant amount of work
    chroms <- as.character(chroms)
    not.in <- chroms[!(chroms %in% seqnames(object@genome.object))]
    if (length(not.in) > 0) {
        cat("the following chromosomes are not recognized:\n")
        print(not.in)
        cat("valid chromosomes names for genome", object@genome.string, 'are:\n')
        print(seqnames(object@genome.object))
        stop('invalid chromosomes, see above for details')
    }

    chrom.refine.records <- setNames(lapply(chroms, abmodel.fit.one.chrom,
        path=path, sc.sample=object@single.cell,
        genome.object=object@genome.object,
        hsnp.tilesize=hsnp.tilesize, n.tiles=n.tiles,
        refine.n.steps=refine.n.steps, n.chunks=n.chunks,
        n.logp.samples.per.chunk=samples.per.chunk), chroms)

    chrom.refine.records
})



# Actually calculate the AB estimates.
# Currently, computes AB at ALL SITES rather than just somatic candidates.
# This will be useful in the future for mosaics and perhaps for using more
# germline hSNPs to model what somatic mutations should look like.
#
# IMPORTANT: AB estimation benefits greatly from chunking. However,
# AB estimation uses a window of 100kb up and downstream from every site at
# which AB is being estimated. For sites at the edge of each chunk, data
# either upstream or downstream will not be available. To solve this,
# training data must be read in AGAIN with the 100kb flanking regions added.
setGeneric("compute.ab.estimates", function(object, n.cores=1, quiet=FALSE)
    standardGeneric("compute.ab.estimates"))
setMethod("compute.ab.estimates", "SCAN2", function(object, n.cores=1, quiet=FALSE) {
    check.slots(object, c('gatk', 'ab.fits'))

    flank <- 1e5 # currently not configurable by user, partly by design

    if (n.cores != 1)
        stop('the n.cores argument is currently unsupported')

    sites <- object@gatk[,.(chr,pos,refnt,altnt)]

    # need to extend region only if this is a chunked object. otherwise,
    # we already have the full table.
    if (!is.null(object@region)) {
        path <- object@integrated.table.path
        if (!quiet) cat('Importing extended hSNP training data from', path, 'using extended range\n')
        extended.range <- GRanges(seqnames=seqnames(object@region)[1],
            ranges=IRanges(start=start(object@region)-flank, end=end(object@region)+flank))
        extended.training.hsnps <- read.training.hsnps(path, sample.id=object@single.cell, region=extended.range, quiet=quiet)
        if (!quiet)
            cat(sprintf("hSNP training sites: %d, extended training sites: %d\n",
                nrow(object@gatk[training.site==TRUE]), nrow(extended.training.hsnps)))
        training.hsnps <- extended.training.hsnps
    } else {
        training.hsnps <- object@gatk[training.site == TRUE]
    }

    # Splitting by chromosome is not for parallelization; the AB model
    # is fit separately for each chromosome and thus applies different
    # parameter estimates.
    chroms <- unique(sites$chr)
    do.work <- function(chrom, sites, ab.fit, hsnps) {
        if (!quiet)
            cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
                nrow(sites), chrom, min(sites$pos), max(sites$pos)))
        time.elapsed <- system.time(z <- infer.gp1(ssnvs=sites, fit=ab.fit,
            hsnps=hsnps, flank=1e5, verbose=!quiet))
        if (!quiet) print(time.elapsed)
        z
    }
    if (length(chroms) == 1) {
        # slightly more efficient for real use cases with chunked computation
        ab <- do.work(chrom=chroms, sites=sites,
            ab.fit=object@ab.fits[chroms,,drop=FALSE],
            hsnps=training.hsnps)
    } else {
        ab <- do.call(rbind, lapply(chroms, function(chrom) {
            hsnps <- training.hsnps[chr == chrom]
            ab.fit <- object@ab.fits[chrom,,drop=FALSE]
            do.work(chrom=chrom, sites=sites, ab.fit=ab.fit, hsnps=hsnps)
        }))
    }

    if (!is.null(ab)) {
        object@gatk[, c('ab', 'gp.mu', 'gp.sd') := 
            list(1/(1+exp(-ab[,'gp.mu'])), ab[,'gp.mu'], ab[,'gp.sd'])]
        object@ab.estimates <- data.frame(sites=nrow(ab))
    } else {
        # Chunks are sometimes empty
        # Add dummy columns so rbind() works with other non-empty chunks
        object@gatk[, c('ab', 'gp.mu', 'gp.sd') := 
            list(numeric(0), numeric(0), numeric(0))]
        object@ab.estimates <- data.frame(sites=0)
    }
    object
})


setGeneric("compute.models", function(object, verbose=TRUE)
    standardGeneric("compute.models"))
setMethod("compute.models", "SCAN2", function(object, verbose=TRUE) {
    check.slots(object, c('gatk', 'ab.estimates'))

    if (nrow(object@gatk) > 0) {
        matched.gp.mu <- match.ab(af=object@gatk$af, gp.mu=object@gatk$gp.mu)
        pvb <- compute.pvs.and.betas(object@gatk$scalt, object@gatk$dp,
                                    matched.gp.mu, object@gatk$gp.sd, verbose=verbose)
    } else {
        pvb <- data.table(abc.pv=numeric(0), lysis.pv=numeric(0), lysis.beta=numeric(0), mda.pv=numeric(0), mda.beta=numeric(0))
    }

    object@gatk[, c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta') := pvb]
    object@mut.models <- data.frame(sites=nrow(pvb))
    object
})


setGeneric("compute.fdr.prior.data", function(object, mode='legacy', quiet=FALSE)
    standardGeneric("compute.fdr.prior.data"))
setMethod("compute.fdr.prior.data", "SCAN2", function(object, mode='legacy', quiet=FALSE) {
    check.slots(object, c('gatk', 'static.filter.params'))

    # FIXME: after integrated table, should probably just use somatic.candidate
    # column here.
    if (mode == 'legacy') {
        # in legacy mode, only candidate sites passing a small set of pre-genotyping
        # crtieria were used.
        bulk.sample <- object@bulk
        bulk.gt <- object@gatk[[bulk.sample]]
        min.sc.alt <- object@static.filter.params$min.sc.alt
        min.sc.dp <- object@static.filter.params$min.sc.dp
        min.bulk.dp <- object@static.filter.params$min.bulk.dp
        max.bulk.alt <- object@static.filter.params$max.bulk.alt
        cand <- object@gatk[
            balt <= max.bulk.alt &
            bulk.gt == '0/0' &
            dbsnp == '.' &
            scalt >= min.sc.alt &
            dp >= min.sc.dp &
            bulk.dp >= min.bulk.dp &
            (is.na(balt.lowmq) | balt.lowmq <= max.bulk.alt)]
        hsnps=object@gatk[training.site == TRUE & scalt >= min.sc.alt]
    } else {
        stop("only the legacy mode is currently implemented. check back soon.")
        # non-legacy mode will apply static filter params, which is almost what is
        # done above.
        if (!('static.filter' %in% colnames(object@gatk)))
            stop('must compute static filters before running compute.fdr.priors')
    }

    # Returns a list of FDR prior data. Also record the mode used.
    object@fdr.prior.data <-
        c(compute.fdr.prior.data.for.candidates(candidates=cand, hsnps=hsnps, random.seed=0, quiet=quiet), mode=mode)
    object
})


setGeneric("compute.fdr", function(object, path, mode='legacy')
    standardGeneric("compute.fdr"))
setMethod("compute.fdr", "SCAN2", function(object, path, mode='legacy') {
    check.slots(object, c('gatk', 'ab.estimates', 'mut.models'))

    if (!missing(path) & !is.null(slot(object, 'fdr.prior.data')))
        stop('fdr.prior.data is already loaded; path to new fdr.prior.data would overwrite')

    if (!missing(path))
        object@fdr.prior.data <- get(load(path))

    # First use NT/NA tables to assign NT and NA to every site
    check.slots(object, 'fdr.prior.data')
    nt.na <- estimate.fdr.priors(object@gatk, object@fdr.prior.data)
    nt <- nt.na$nt
    na <- nt.na$na
    object@gatk[, c('nt', 'na') := list(nt, na)]

    # Next compute lysis.fdr and mda.fdr, which represent the false discovery rate
    # of a population of candidate mutation sites with the same DP and VAF as the
    # site in question.
    # Legacy computation (finding min FDR over all alphas) and legacy candidate set.
    # Legacy computation is too slow for applying to all sites.
    if (mode == 'legacy') {
        bulk.sample <- object@bulk
        bulk.gt <- object@gatk[[bulk.sample]]
        cand <- object@gatk[
            balt == 0 &
            bulk.gt == '0/0' &
            dbsnp == '.' &
            scalt >= object@static.filter.params$min.sc.alt &
            dp >= object@static.filter.params$min.sc.dp]
            # legacy did NOT require passing min.bulk.dp or 0 bulk alt reads at low MQ
            #bulk.dp >= object@static.filter.params$min.bulk.dp]
            #(is.na(balt.lowmq) | balt.lowmq == 0)]

        matched.gp.mu <- match.ab(af=cand$af, gp.mu=cand$gp.mu)
        cand[, c('lysis.fdr', 'mda.fdr') :=
            compute.fdr.legacy(altreads=scalt, dp=dp,
                gp.mu=matched.gp.mu, gp.sd=gp.sd, nt=nt, na=na)]
        object@gatk[cand, on=.(chr,pos,refnt,altnt),
            c('lysis.fdr', 'mda.fdr') := list(i.lysis.fdr, i.mda.fdr)]
        object@fdr <- list(mode=mode, sites=nrow(cand))
    } else if (mode == 'new') {
        # XXX: TODO: calculate adjusted NA/NT values for hSNPs. Equivalent to
        # previous leave-one-out approaches, because AB estimation always leaves
        # out the site being estimated (whether hSNP or somatic candidate).
        #
        # HOWEVER, this FDR heuristic on hSNPs IS NOT a good way to actually call
        # germline hSNPs if that is your goal. These FDR heuristics are tuned to
        # the specific set of CANDIDATE SOMATIC MUTATIONS detected for this cell.
        object@gatk[, c('lysis.fdr', 'mda.fdr') :=
            list(lysis.pv*na / (lysis.pv*na + lysis.beta*nt),
                 mda.pv*na / (mda.pv*na + mda.beta*nt))]
        object@fdr <- list(mode=mode, sites=nrow(object@gatk))
    } else
        stop(sprintf("unrecognized mode '%s', expecting either 'legacy' or 'new'", mode))

    object
})


cigar.emp.score <- function (training, test, which = c("id", "hs"), quiet=FALSE) {
    xt <- training[[paste0(which, ".score.x")]]
    yt <- training[[paste0(which, ".score.y")]]
    x <- test[[paste0(which, ".score.x")]]
    y <- test[[paste0(which, ".score.y")]]

    dp <- test$dp.cigars
    bulkdp <- test$dp.cigars.bulk

    progressr::with_progress({
        if (!quiet) p <- progressr::progressor(along=1:(length(x)/100))
        ret <- future.apply::future_mapply(function(dp, bulkdp, x, y, i) {
            if (!quiet & i %% 100 == 1) p()
            ifelse(dp == 0 | bulkdp == 0, 0,
                mean(xt >= x & yt >= y, na.rm = T))
        }, test$dp.cigars, test$dp.cigars.bulk, x, y, 1:length(x))
    }, enable=!quiet)

    # future_mapply returns list() when x is length 0. make this a 0-length
    # numeric so quantile() doesn't fail on it later.
    if (length(ret) == 0)
        return(numeric(0))
    else
        return(ret)
}


compute.cigar.scores <- function(cigar.data) {
    cigar.data[, c('id.score.y', 'id.score.x', 'hs.score.y', 'hs.score.x') :=
        list(ID.cigars / dp.cigars,
             ID.cigars.bulk / dp.cigars.bulk,
             HS.cigars / dp.cigars,
             HS.cigars.bulk / dp.cigars.bulk)]
}


# modifies 'data' by reference
compute.excess.cigar <- function(data, cigar.training, quiet=FALSE) {
    pc <- perfcheck("excess indel (I/D) ops",
            idopscores <- cigar.emp.score(training=cigar.training, test=data, which='id', quiet=quiet),
        report.mem=FALSE),
    if (!quiet) cat(pc, '\n')
    pc <- perfcheck("excess clipping (H/S) ops",
            hsopscores <- cigar.emp.score(training=cigar.training, test=data, which='hs', quiet=quiet),
        report.mem=FALSE)
    if (!quiet) cat(pc, '\n')

    data[, c('id.score', 'hs.score') := list(idopscores, hsopscores)]
}


read.cigar.data <- function(path, region, quiet=FALSE) {
    if (!quiet) cat('Importing CIGAR stats from', path, '\n')
    col.classes <- c('character', rep('integer', 6))
    read.tabix.data(path=path, region=region, quiet=quiet, colClasses=col.classes)
}

        
setGeneric("add.cigar.data", function(object, sc.cigars.path, bulk.cigars.path, quiet=FALSE)
    standardGeneric("add.cigar.data"))
setMethod("add.cigar.data", "SCAN2", function(object, sc.cigars.path, bulk.cigars.path, quiet=FALSE) {
    check.slots(object, 'gatk')

    sc <- read.cigar.data(sc.cigars.path, region=object@region, quiet=quiet)
    bulk <- read.cigar.data(bulk.cigars.path, region=object@region, quiet=quiet)

    if (!quiet) cat('joining CIGAR data..\n')
    object@gatk[sc, on=c('chr', 'pos'),
        c('M.cigars', 'ID.cigars', 'HS.cigars', 'other.cigars', 'dp.cigars') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]
    object@gatk[bulk, on=c('chr', 'pos'),
        c('M.cigars.bulk', 'ID.cigars.bulk', 'HS.cigars.bulk', 'other.cigars.bulk', 'dp.cigars.bulk') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]

    if (!quiet) cat('computing CIGAR op rates..\n')
    compute.cigar.scores(object@gatk)  # modifies by reference
    object@cigar.data <- data.frame(sc.sites=nrow(sc), bulk.sites=nrow(bulk),
        sc.path=sc.cigars.path, bulk.path=bulk.cigars.path)
    object
})


cigar.get.null.sites <- function(object, path=NULL, legacy=TRUE, quiet=FALSE) {
    if (is.null(path)) {
        check.chunked(object,
            'excess cigar scores should be computed on full chr1-chr22 data, not chunked data')

        if (legacy) {
            check.slots(object, c('gatk', 'cigar.data'))
            null.sites <- object@gatk[resampled.training.site==TRUE]
            if (!quiet) cat(sprintf('LEGACY: computing CIGAR op rates only at resampled training sites (n=%d)..\n',
                nrow(null.sites)))
        } else {
            check.slots(object, c('gatk', 'cigar.data'))
            null.sites <- object@gatk[training.site==TRUE]
            if (!quiet) {
                cat(sprintf('computing CIGAR op rates for all training sites (n=%d)..\n',
                    nrow(null.sites)))
                cat('WARNING: using the full set of hSNPs may be prohibitively slow\n')
            }
        }
    } else {
        null.sites <- data.table::fread(path)
    }

    null.sites
}


setGeneric("compute.excess.cigar.scores", function(object, path=NULL, legacy=TRUE, quiet=FALSE)
    standardGeneric("compute.excess.cigar.scores"))
setMethod("compute.excess.cigar.scores", "SCAN2", function(object, path=NULL, legacy=TRUE, quiet=FALSE) {
    null.sites <- cigar.get.null.sites(object, path, legacy, quiet)
    compute.excess.cigar(data=object@gatk, cigar.training=null.sites, quiet=quiet)
    object@excess.cigar.scores <- data.frame(training.sites=nrow(null.sites), legacy=legacy)
    object
})


setGeneric("add.static.filter.params", function(object, config.path,
    min.sc.alt=2, min.sc.dp=6,
    max.bulk.alt=0, min.bulk.dp=11, exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
        standardGeneric("add.static.filter.params"))
setMethod("add.static.filter.params", "SCAN2",
function(object, config.path,
    min.sc.alt=2, min.sc.dp=6, max.bulk.alt=0, min.bulk.dp=11,
    exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
{
    if (!missing(config.path)) {
        yaml <- yaml::read_yaml(config.path)
        min.sc.alt <- yaml$min_sc_alt
        min.sc.dp <- yaml$min_sc_dp
        max.bulk.alt <- yaml$max_bulk_alt
        min.bulk.dp <- yaml$min_bulk_dp
        # exclude.dbsnp, cg.id.q and cg.hs.q are not user configurable at the moment
    }
    object@static.filter.params <- list(
        min.sc.alt=min.sc.alt, min.sc.dp=min.sc.dp,
        max.bulk.alt=max.bulk.alt, min.bulk.dp=min.bulk.dp,
        exclude.dbsnp=exclude.dbsnp, cg.id.q=cg.id.q, cg.hs.q=cg.hs.q)
    object
})
        

setGeneric("compute.static.filters", function(object, exclude.dbsnp=TRUE)
        standardGeneric("compute.static.filters"))
setMethod("compute.static.filters", "SCAN2", function(object, exclude.dbsnp=TRUE) {
    check.slots(object, c('gatk', 'cigar.data', 'mut.models'))

    qid <- quantile(object@gatk[training.site == TRUE]$id.score,
        prob=object@static.filter.params$cg.id.q, na.rm=TRUE)
    qhs <- quantile(object@gatk[training.site == TRUE]$hs.score,
        prob=object@static.filter.params$cg.hs.q, na.rm=TRUE)

    # having some issues with data.table accessing objects/dataframes inside the expression
    max.bulk.alt <- object@static.filter.params$max.bulk.alt
    min.sc.dp <- object@static.filter.params$min.sc.dp
    min.bulk.dp <- object@static.filter.params$min.bulk.dp
    min.sc.alt <- object@static.filter.params$min.sc.alt
    object@gatk[, c('cigar.id.test', 'cigar.hs.test', 'lowmq.test',
            'dp.test', 'abc.test', 'min.sc.alt.test', 'max.bulk.alt.test',
            'dbsnp.test') :=
                list(id.score > qid,
                     hs.score > qhs,
                     is.na(balt.lowmq) | balt.lowmq <= max.bulk.alt,
                     dp >= min.sc.dp & bulk.dp >= min.bulk.dp,
                     abc.pv > 0.05,
                     scalt >= min.sc.alt,
                     balt <= max.bulk.alt,
                     !exclude.dbsnp | dbsnp == '.')]
    object@gatk[,static.filter :=
        cigar.id.test & cigar.hs.test & lowmq.test & dp.test &
        abc.test & min.sc.alt.test & max.bulk.alt.test & dbsnp.test]
    object
})


# When reading the integrated table: read all metadata columns (currently 1-18)
# and only the 3 genotype/count columns corresponding to `sample.id`. This can
# save significant memory overhead in large (100+ cell) projects.
#
# Next, annotate the integrated table with information corresponding to
# `sample.id`. This includes training site definitions (which varies from cell
# to cell depending on whether there is no data at the site in that single cell;
# gt=./.) and assigning single cell ref/alt read counts to phased haplotypes.
#
# N.B. parsimony phasing (adjust.phase()) used to be called here, but that is not
# a good idea. That should happen in the make.integrated.table pipeline where the
# full single cell+bulk count table is available. Information shared across single
# cells is useful to improve phasing and the final phase decision should be
# consistent across single cells.
read.and.annotate.integrated.table <- function(path, sample.id, region=NULL, quiet=FALSE) {
    tr <- read.table.1sample(path, sample.id, n.meta.cols=18, region=region, quiet=quiet)
    setindex(tr, resampled.training.site)

    # Add some convenient calculations
    tr[, dp := scalt + scref]
    tr[, af := scalt / dp]
    data.table::setkey(tr, chr, pos, refnt, altnt)
    setindex(tr, muttype) # allow for fast selection of SNVs or indels

    sc.gt <- tr[[sample.id]]  # the column named after the sample is the GATK GT string for that sample
    tr[, training.site := (phased.gt == '1|0' | phased.gt == '0|1') & sc.gt != './.' & bulk.gt != './.']
    tr[, c('phased.hap1', 'phased.hap2') :=
        list(ifelse(phased.gt == '0|1', scref, scalt),
             ifelse(phased.gt == '0|1', scalt, scref))]
    tr
}


read.training.hsnps <- function(path, sample.id, region=NULL, quiet=FALSE) {
    read.and.annotate.integrated.table(path=path, sample.id=sample.id, region=region, quiet=quiet)[training.site == TRUE & muttype == 'snv']
}
