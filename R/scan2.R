setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClassUnion('null.or.BSgenome', c('NULL', 'BSgenome'))
setClassUnion('null.or.GRanges', c('NULL', 'GRanges'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.logical', c('NULL', 'logical'))

# if adding new slots to the class, be sure to update concat(...)
# to combine the slots properly.
setClass("SCAN2", slots=c(
    region='null.or.GRanges',
    genome.string='character',
    genome.object='null.or.BSgenome',
    single.cell='character',
    bulk='character',
    gatk="null.or.df",
    training.data="null.or.df",
    resampled.training.data="null.or.list",
    ab.fits='null.or.df',
    static.filter.params='null.or.list',
    gatk.lowmq="null.or.df",
    ab.estimates='null.or.df',
    mut.models='null.or.df',
    cigar.data='null.or.df',
    excess.cigar.scores='null.or.df',
    static.filters='null.or.df',
    static.filter='null.or.logical',
    fdr.priors='null.or.df',
    fdrs='null.or.df'))


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
    new("SCAN2", single.cell=single.cell, bulk=bulk,
        genome.string=genome,
        genome.object=genome.string.to.bsgenome.object(genome),
        region=region,
        gatk=NULL,
        training.data=NULL,
        ab.fits=NULL,
        # these slots used to hold tables; now their data is incorporated into @gatk and
        # they only record analysis parameters, if any apply.
        gatk.lowmq=NULL,
        ab.estimates=NULL,
        mut.models=NULL,
        cigar.data=NULL,
        excess.cigar.scores=NULL,
        static.filter.params=NULL,
        fdr.priors=NULL,
        fdrs=NULL)
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
        if (length(object@region) != 1) {
            stop('"region" must contain exactly 1 range')
        }
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

    if (!is.null(object@gatk.lowmq)) {
        if (!all(c('lowmq.scref', 'lowmq.scalt', 'lowmq.bref', 'lowmq.balt') %in%
            colnames(object@gatk)))
            return("low mapping quality GATK output is missing from @gatk")
    }

    if (!is.null(object@training.data)) {
        if (!all(c('hap1', 'hap2', 'phgt', 'training.site.full') %in%
            colnames(object@gatk)))
            return("training.data was not properly imported")
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
            if (s == 'gatk.lowmq')
                cat("must import low mapping quality GATK read counts first (see: read.gatk.lowmq())\n")
            if (s == 'training.data')
                cat("must import hSNP training data first (see: add.training.data())\n")
            if (s == 'resampled.training.data')
                cat("must resample hSNP training data first (see: resample.training.data())\n")
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
            if (s == 'fdr.priors')
                cat("must compute FDR priors first (see: compute.fdr.priors())\n")
        }
    }
    if (error.occurred & abort)
        stop("One or more required slots are missing. See above for details.")
    return(error.occurred)
}


setMethod("show", "SCAN2", function(object) {
    cat("#", is(object)[[1]], "\n")
    if (!is.null(object@region)) {
        cat("#   Region:",
            paste0(as.character(seqnames(object@region)[1]), ':',
            start(object@region)[1], '-',
            end(object@region)[1]), '\n')
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

    cat("#   GATK low mapping quality:",
        ifelse(is.null(object@gatk.lowmq), " (no data)", paste(object@gatk.lowmq$sites, 'sites')), '\n')

    cat("#   AB model training hSNPs:")
    if (is.null(object@training.data)) {
        cat(" (no data)\n")
    } else {
        cat('', object@training.data$sites,
            sprintf("phased sites (hap1=%d, hap2=%d)",
                sum(object@gatk$training.phgt=='1|0', na.rm=TRUE),
                sum(object@gatk$training.phgt=='0|1', na.rm=TRUE)),"\n")
        zzz <- do.call(rbind, lapply(split(object@gatk, object@gatk$chr),
            function(td) {
                td$af <- td$training.hap1/(td$training.hap1+td$training.hap2)
                cbind(dist=diff(td$pos), af1=td$af[-nrow(td)], af2=td$af[-1])
        }))
        cors <- sapply(10^(2:5), function(threshold) {
            z <- zzz[zzz[,1] <= threshold,]; cor(z[,2], z[,3], use='complete.obs') }) 
        cat('#       VAF correlation between neighboring hSNPs:\n')
        cat('#           <100 bp', round(cors[1], 3),
            '<1000 bp', round(cors[2], 3),
            '<10 kbp', round(cors[3], 3),
            '<100 kbp', round(cors[4], 3), '\n')
        if (!is.null(object@resampled.training.data)) {
            cat('#        ', nrow(object@gatk[resampled.training.site == TRUE]),
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
        if (!is.null(object@training.data)) {
            xs <- round(object@gatk[training.site==TRUE,
                .(mean=mean(gp.mu), cor=cor(af, ab, use='complete.obs'))],3)
            cat('#       mean at training hSNPs:', xs$mean, '\n')
            cat('#       correlation with VAF at training hSNPs', xs$cor, '\n')
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
    if (is.null(object@static.filters)) {
        cat("(not applied)\n")
    } else {
        cat(sum(object@static.filter, na.rm=TRUE), "retained",
            sum(!object@static.filter, na.rm=TRUE), "removed",
            sum(is.na(object@static.filter)), "NA\n")
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

    # XXX: TODO: would be nice a nice sanity check to make sure each SCAN2
    # class in the list is compatible (i.e., from the same reference genome/
    # sample IDs)
    init <- args[[1]]

    ret <- make.scan(init@single.cell, init@bulk, genome=init@genome.string)
    # rbindlist is a special data.table function for quickly 
    ret@gatk <- rbindlist(lapply(args, function(a) a@gatk))

    ret@training.data <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@training.data), 0, a@training.data$sites))))

    if (any(sapply(args, function(a) !is.null(a@resampled.training.data)))) {
        warning('concat discards information used to resample training data; resampled sites are, however, retained')
        ret@resampled.training.data <- list(sites=sum(sapply(args, function(a) ifelse(is.null(a@resampled.training.data), 0, a@resampled.training.data$selection$keep))))
    } else {
        ret@resampled.training.data <- NULL
    }

    # policy: ab.fits has to be the same for all chunks being concat()ed
    # if that's true, just use the first one
    if (any(sapply(args, function(a) any(as.matrix(init@ab.fits) != as.matrix(a@ab.fits)))))
        stop('@ab.fits must be identical for all concat() elements')
    ret@ab.fits <- init@ab.fits

    # policy: same as above: all static filter params must be identical
    if (any(sapply(args, function(a) any(as.vector(init@static.filter.params) != as.vector(a@static.filter.params)))))
        stop('@static.filter.params must be identical for all concat() elements')
    ret@static.filter.params <- init@static.filter.params

    ret@gatk.lowmq <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@gatk.lowmq), 0, a@gatk.lowmq$sites))))

    ret@ab.estimates <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@ab.estimates), 0, a@ab.estimates$sites))))

    ret@mut.models <- data.frame(sites=sum(sapply(args, function(a) ifelse(is.null(a@mut.models), 0, a@mut.models$sites))))

    ret@cigar.data <- data.frame(sc.sites=sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$sc.sites))),
                               bulk.sites=sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$bulk.sites)))
    )

    ret@excess.cigar.scores <- data.frame(training.sites=sum(sapply(args, function(a) ifelse(is.null(a@excess.cigar.scores), 0, a@excess.cigar.scores$training.sites))))

    ret@fdr.priors <- NULL
    ret@fdrs <- NULL

    ret
})



# Some extra work to make sure we only read in the part of the
# table relevant to these two samples. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
read.gatk.table.2sample <- function(path, sc.sample, bulk.sample, region) {
    cat("Importing GATK table..\n")

    # Step 1: just get the header and detect the columns corresponding to
    # sc.sample and bulk.sample to avoid reading the full matrix.
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- sub('^#', '', Rsamtools::headerTabix(tf)$header) # strip the leading #
    col.strings <- strsplit(header, '\t')[[1]]
    tot.cols <- length(col.strings)
    sc.idx <- which(col.strings == sc.sample)
    bulk.idx <- which(col.strings == bulk.sample)
    cat("Selecting columns:\n")
    for (i in 1:length(col.strings)) {
        if (i <= 7) {
            cat(sprintf("    (%d)", i), col.strings[i], '\n')
        } else if (any(i == sc.idx + 0:2)) {
            cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
        } else if (any(i == bulk.idx + 0:2)) {
            cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
        }
    }
    
    # Step 2: really read the tables in, but only the relevant columns
    cols.to.read <- rep("NULL", tot.cols)
    # First 7 are chr, pos, dbsnp, refnt, altnt, mq, mqrs
    cols.to.read[1:7] <- c('character', 'integer', rep('character', 3), 'numeric', 'numeric')
    # Read 3 columns for the single cell, 3 columns for bulk
    cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
    cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')
    if (is.null(region))
        gatk <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf)[[1]]), colClasses=cols.to.read)
    else
        gatk <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf, param=region)[[1]]), colClasses=cols.to.read)
    cat("Read", nrow(gatk), 'lines\n')
    close(tf)
    new.sc.idx <- which(colnames(gatk) == sc.sample)
    new.bulk.idx <- which(colnames(gatk) == bulk.sample)
    colnames(gatk)[new.sc.idx+1:2] <- c('scref', 'scalt')
    colnames(gatk)[new.bulk.idx+1:2] <- c('bref', 'balt')

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- c(col.strings[1:7], sc.sample, c('scref', 'scalt'), bulk.sample, c('bref', 'balt'))
    gatk <- gatk[,..cols.to.keep]

    gatk
}


setGeneric("read.gatk", function(object, path)
    standardGeneric("read.gatk"))
setMethod("read.gatk", "SCAN2", function(object, path) {
    gatk <- read.gatk.table.2sample(path, object@single.cell, object@bulk, object@region)

    # Add some convenient calculations
    gatk$dp <- gatk$scalt + gatk$scref
    gatk$af <- gatk$scalt / gatk$dp
    gatk$bulk.dp <- gatk$balt + gatk$bref
    gatk$bulk.af <- gatk$balt / gatk$bulk.dp
    data.table::setkey(gatk, chr, pos, refnt, altnt)

    # Determine SNV/indel status and then annotate mutation signature channels
    gatk[, muttype := ifelse(nchar(refnt) == 1 & nchar(altnt) == 1, 'snv', 'indel')]
    # allow for fast selection of SNVs or indels
    setindex(gatk, muttype)

    gatk[muttype == 'snv',
        mutsig := get.3mer(chr=chr, pos=pos, refnt=refnt, altnt=altnt, genome=object@genome.object)]
    chs <- classify.indels(gatk[muttype == 'indel'], genome.string=object@genome.string)
    gatk[muttype == 'indel', mutsig := chs]

    object@gatk <- gatk
    object
})


setGeneric("read.gatk.lowmq", function(object, path)
    standardGeneric("read.gatk.lowmq"))
setMethod("read.gatk.lowmq", "SCAN2", function(object, path) {
    check.slots(object, 'gatk')

    lowmq <- read.gatk.table.2sample(path, object@single.cell, object@bulk, object@region)
    data.table::setkey(lowmq, chr, pos, refnt, altnt)

    object@gatk[lowmq, on=.(chr,pos,refnt,altnt), c('scref.lowmq', 'scalt.lowmq', 'bref.lowmq', 'balt.lowmq') := list(i.scref, i.scalt, i.bref, i.balt)]
    object@gatk.lowmq <- data.frame(sites=nrow(lowmq), path=path)
    object
})


# Add AB Gaussian process parameter fits for each chromosome.
# These take a long time to compute, so the snakemake pipeline will almost
# always be necessary here.
setGeneric("add.ab.fits", function(object, path)
    standardGeneric("add.ab.fits"))
setMethod("add.ab.fits", "SCAN2", function(object, path) {
    fitlist <- get(load(path))
    object@ab.fits <- do.call(rbind, fitlist)
    object
})



# Actually calculate the AB estimates.
# Currently, computes AB at ALL SITES.
setGeneric("compute.ab.estimates", function(object, n.cores=1)
    standardGeneric("compute.ab.estimates"))
setMethod("compute.ab.estimates", "SCAN2", function(object, n.cores=1) {
    check.slots(object, c('gatk', 'training.data', 'ab.fits'))

    if (n.cores != 1)
        stop('the n.cores argument is currently unsupported')

    sites <- object@gatk[,.(chr,pos,refnt,altnt)]

    # Splitting by chromosome is not for parallelization; the AB model
    # is fit separately for each chromosome and thus applies different
    # parameter estimates.
    chroms <- unique(sites$chr)
    do.work <- function(chrom, sites, ab.fit, hsnps) {
        print(ab.fit)
        cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
            nrow(sites), chrom, min(sites$pos), max(sites$pos)))
cat("---------------- FIXME FIXME FIXME --------------------\n")
cat("---------------- GENERATING RANDOM AB ESTS --------------------\n")
cat("---------------- BECAUSE LAPACKE NOT INSTALLED--------------------\n")
        return(data.frame(gp.mu=rnorm(nrow(sites)), gp.sd=rgamma(nrow(sites), shape=1)))
        time.elapsed <- system.time(z <- infer.gp1(ssnvs=sites, fit=ab.fit,
            hsnps=hsnps, flank=1e5, verbose=TRUE))
        print(time.elapsed)
        z
    }
    if (length(chroms) == 1) {
        # slightly more efficient for real use cases with chunked computation
        ab <- do.work(chrom=chroms, sites=sites,
            ab.fit=object@ab.fits[chroms,,drop=FALSE],
            hsnps=object@gatk[training.site == TRUE,])
    } else {
        ab <- do.call(rbind, lapply(chroms, function(chrom) {
            hsnps <- object@gatk[chr == chrom & training.site == TRUE,]
            ab.fit <- object@ab.fits[chrom,,drop=FALSE]
            do.work(chrom=chrom, sites=sites, ab.fit=ab.fit, hsnps=hsnps)
        }))
    }

    object@gatk[, c('ab', 'gp.mu', 'gp.sd') := 
        list(1/(1+exp(-ab$gp.mu)), ab$gp.mu, ab$gp.sd)]
    object@ab.estimates <- data.frame(sites=nrow(ab))
    object
})


setGeneric("compute.models", function(object)
    standardGeneric("compute.models"))
setMethod("compute.models", "SCAN2", function(object) {
    check.slots(object, c('gatk', 'ab.estimates'))

    pvb <- compute.pvs.and.betas(object@gatk$scalt, object@gatk$dp,
                                 object@gatk$gp.mu, object@gatk$gp.sd)
    object@gatk[, c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta') := pvb]
    object@mut.models <- data.frame(sites=nrow(pvb))
    object
})


setGeneric("compute.fdr.priors", function(object, mode='legacy')
    standardGeneric("compute.fdr.priors"))
setMethod("compute.fdr.priors", "SCAN2", function(object, mode) {
    check.slots(object, c('gatk', 'static.filter'))
    # minimized objects have training.data deleted, but training.site is
    # annotated in the gatk data frame.
    if (!('training.site' %in% colnames(object@gatk)))
        check.slots(object, 'training.data')

    cat("*** FIXME: filtering out static.filter=NA; should not contain NAs.\n")
    if (mode == 'legacy') {
        # in legacy mode, all candidate sites passing a small set of pre-genotyping
        # crtieria were used.
        cand <- object@gatk[
            object@gatk$balt == 0 &
            object@gatk[,11] == '0/0' &
            object@gatk$dbsnp == '.' &
            object@gatk$scalt >= object@static.filter.params$min.sc.alt &
            object@gatk$dp >= object@static.filter.params$min.sc.dp &
            object@gatk$bulk.dp >= object@static.filter.params$min.bulk.dp &
            (is.na(object@gatk.lowmq$balt) | object@gatk.lowmq$balt == 0),]
        hsnps=object@gatk[
            object@gatk$training.site &
            object@gatk$scalt >= object@static.filter.params$min.sc.alt,]
    } else {
        stop("only the legacy mode is currently implemented. check back soon.")
    }

    cat(nrow(cand), 'somatic candidates\n')
    cat(nrow(hsnps), 'hsnps\n')
    fdr.priors <-
        estimate.fdr.priors(candidates=cand, hsnps=hsnps, random.seed=0)
    rownames(fdr.priors) <- rownames(cand)

    join.cols <- c('chr','pos','refnt','altnt')
    object@fdr.priors <- plyr::join(object@gatk[,join.cols],
        cbind(cand[,join.cols], fdr.priors))[,-(1:4)]
    object
})


setGeneric("compute.fdr", function(object, mode='legacy')
    standardGeneric("compute.fdr"))
setMethod("compute.fdr", "SCAN2", function(object, mode='legacy') {
    check.slots(object, c('gatk', 'ab.estimates', 'mut.models', 'fdr.priors'))

    # Legacy computation (finding min FDR over all alphas) and legacy candidate set.
    # Importantly, the legacy method is far too slow to apply to all genomic loci,
    # so extra work must be done to only compute over the legacy set.
    if (mode == 'legacy') {
        cand <-
            object@gatk$balt == 0 &
            object@gatk[,11] == '0/0' &
            object@gatk$dbsnp == '.' &
            object@gatk$scalt >= object@static.filter.params$min.sc.alt &
            object@gatk$dp >= object@static.filter.params$min.sc.dp &
            object@gatk$bulk.dp >= object@static.filter.params$min.bulk.dp &
            (is.na(object@gatk.lowmq$balt) | object@gatk.lowmq$balt == 0)

        fdrs <- compute.fdr.legacy(
            altreads=object@gatk$scalt[cand],
            dp=object@gatk$dp[cand],
            gp.mu=object@ab.estimates$gp.mu[cand],
            gp.sd=object@ab.estimates$gp.sd[cand],
            nt=object@fdr.priors$nt[cand],
            na=object@fdr.priors$na[cand])

        join.cols <- c('chr','pos','refnt','altnt')
        object@fdrs <- plyr::join(object@gatk[,join.cols],
            cbind(object@gatk[cand,join.cols], fdrs))[,-(1:4)]
    } else if (mode == 'new')
        object@fdrs <- compute.fdr.new(object@mut.models, object@fdr.priors)
    else
        stop(sprintf("unrecognized mode '%s', expecting either 'legacy' or 'new'", mode))

    object
})


cigar.emp.score <- function (training, test, which = c("id", "hs")) {
    xt <- training[[paste0(which, ".score.x")]]
    yt <- training[[paste0(which, ".score.y")]]
    x <- test[[paste0(which, ".score.x")]]
    y <- test[[paste0(which, ".score.y")]]

    pbapply::pbmapply(function(xi, yi, bulkdp, dp)
        ifelse(dp == 0 | bulkdp == 0, 0,
            mean(xt >= xi & yt >= yi, na.rm = T)),
        x, y, test$dp.cigars.bulk, test$dp.cigars)
}


compute.cigar.scores <- function(cigar.data) {
    cigar.data[, c('id.score.y', 'id.score.x', 'hs.score.y', 'hs.score.x') :=
        list(ID.cigars / dp.cigars,
             ID.cigars.bulk / dp.cigars.bulk,
             HS.cigars / dp.cigars,
             HS.cigars.bulk / dp.cigars.bulk)]
}


# modifies 'data' by reference
compute.excess.cigar <- function(data, cigar.training) {
    cat("excess indel (I/D) ops..\n")
    data[, c('id.score', 'hs.score') := list(
        cigar.emp.score(training=cigar.training, test=data, which='id'),
        cigar.emp.score(training=cigar.training, test=data, which='hs'))]
}


read.cigar.data <- function(path, region) {
    cat('Importing CIGAR stats from', path, '\n')
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- sub('^#', '', Rsamtools::headerTabix(tf)$header) # strip the leading #
    col.classes <- c('character', rep('integer', 6))
    if (!is.null(region)) {
        tab <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf, param=region)[[1]]),
            colClasses=col.classes)
    } else {
        tab <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf)[[1]]),
            colClasses=col.classes)
    }
    cat("Read", nrow(tab), 'lines\n')
    close(tf)
    tab
}

        
setGeneric("add.cigar.data", function(object, sc.cigars.path, bulk.cigars.path)
    standardGeneric("add.cigar.data"))
setMethod("add.cigar.data", "SCAN2", function(object, sc.cigars.path, bulk.cigars.path) {
    check.slots(object, 'gatk')

    sc <- read.cigar.data(sc.cigars.path, region=object@region)
    bulk <- read.cigar.data(bulk.cigars.path, region=object@region)

    cat('joining CIGAR data to GATK sites..\n')
    object@gatk[sc, on=c('chr', 'pos'),
        c('M.cigars', 'ID.cigars', 'HS.cigars', 'other.cigars', 'dp.cigars') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]
    object@gatk[bulk, on=c('chr', 'pos'),
        c('M.cigars.bulk', 'ID.cigars.bulk', 'HS.cigars.bulk', 'other.cigars.bulk', 'dp.cigars.bulk') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]

    cat('computing CIGAR op rates for GATK sites..\n')
    compute.cigar.scores(object@gatk)  # modifies by reference
    object@cigar.data <- data.frame(sc.sites=nrow(sc), bulk.sites=nrow(bulk),
        sc.path=sc.cigars.path, bulk.path=bulk.cigars.path)
    object
})


setGeneric("compute.excess.cigar.scores", function(object, legacy=TRUE)
    standardGeneric("compute.excess.cigar.scores"))
setMethod("compute.excess.cigar.scores", "SCAN2", function(object, legacy=TRUE) {
    if (legacy) {
        check.slots(object, c('gatk', 'training.data', 'cigar.data', 'resampled.training.data'))
        null.sites <- object@gatk[resampled.training.site==TRUE]
        cat(sprintf('LEGACY: computing CIGAR op rates only at resampled training sites (n=%d)..\n',
            nrow(null.sites)))
    } else {
        check.slots(object, c('gatk', 'training.data', 'cigar.data'))
        null.sites <- object@gatk[training.site==TRUE]
        cat(sprintf('computing CIGAR op rates for all training sites (n=%d)..\n',
            nrow(null.sites)))
        cat('WARNING: using the full set of hSNPs may be prohibitively slow\n')
    }

    compute.excess.cigar(data=object@gatk, cigar.training=null.sites)
    object@excess.cigar.scores <- data.frame(training.sites=nrow(null.sites), legacy=legacy)
    object
})


setGeneric("add.static.filter.parameters", function(object, min.sc.alt=2, min.sc.dp=6,
    max.bulk.alt=0, min.bulk.dp=11, exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
        standardGeneric("add.static.filter.parameters"))
setMethod("add.static.filter.parameters", "SCAN2",
function(object, min.sc.alt=2, min.sc.dp=6, max.bulk.alt=0, min.bulk.dp=11,
    exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
{
    object@static.filter.parameters <- list(
        min.sc.alt=min.sc.alt, min.sc.dp=min.sc.dp,
        max.bulk.alt=max.bulk.alt, min.bulk.dp=min.bulk.dp,
        exclude.dbsnp=exclude.dbsnp, cg.id.q=cg.id.q, cg.hs.q=cg.hs.q)
    object
})
        

setGeneric("compute.static.filters", function(object, exclude.dbsnp=TRUE)
        standardGeneric("compute.static.filters"))
setMethod("compute.static.filters", "SCAN2", function(object, exclude.dbsnp=TRUE) {
    check.slots(object, c('gatk', 'cigar.data', 'gatk.lowmq', 'mut.models'))

    qid <- quantile(object@gatk[training.site == TRUE]$id.score, prob=cg.id.q, na.rm=TRUE)
    qhs <- quantile(object@gatk[training.site == TRUE]$hs.score, prob=cg.hs.q, na.rm=TRUE)

    # having some issues with data.table accessing objects/dataframes inside the expression
    max.bulk.alt <- gatk@static.filter.params$max.bulk.alt
    min.sc.dp <- gatk@static.filter.params$min.sc.dp
    min.bulk.dp <- gatk@static.filter.params$min.bulk.dp
    min.sc.alt <- gatk@static.filter.params$min.sc.alt
    object@gatk[, c('cigar.id.test', 'cigar.hs.test', 'lowmq.test',
            'dp.test', 'abc.test', 'min.sc.alt.test', 'max.bulk.alt.test',
            'dbsnp.test') :=
                list(id.score > qid,
                     hs.score > qhs,
                     is.na(balt.lowmq) | gatk.lowmq$balt <= max.bulk.alt,
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


setGeneric("resample.training.data", function(object, M=20, seed=0) 
        standardGeneric("resample.training.data"))
setMethod("resample.training.data", "SCAN2", function(object, M=20, seed=0) {
    check.slots(object, c('gatk', 'training.data'))

    ret <- resample.hsnps(
        sites=object@gatk[training.site == FALSE & scalt >= 2 & dp >= 6 & bulk.dp >= 11 & balt == 0,],
        hsnps=object@gatk[training.site == TRUE], M=M, seed=seed)

    object@gatk[training.site == TRUE, resampled.training.site := ret$selection$keep]
    object@gatk[is.na(resampled.training.site), resampled.training.site := FALSE]
    setindex(object@gatk, resampled.training.site)
    object@resampled.training.data <- ret
    object
})


read.training.data <- function(path, region=NULL) {
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- sub('^#', '', Rsamtools::headerTabix(tf)$header) # strip the leading #
    col.strings <- strsplit(header, '\t')[[1]]
    # chr, pos, refnt, altnt, gt, hap1, hap2, dp, phgt
    col.classes <- c('character', 'integer', 'character', 'character', 'character', 'integer', 'integer', 'integer', 'character')
    if (is.null(region)) {
        hsnps <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf)[[1]]),
            colClasses=col.classes)
    } else {
        hsnps <- data.table::fread(text=c(header, Rsamtools::scanTabix(tf, param=region)[[1]]),
            colClasses=col.classes)
    }
    close(tf)
    data.table::setkey(hsnps, chr, pos, refnt, altnt)
}


setGeneric("add.training.data", function(object, path)
        standardGeneric("add.training.data"))
setMethod("add.training.data", "SCAN2", function(object, path) {
    cat('Importing hSNP training data from', path, '\n')
    hsnps <- read.training.data(path, object@region)
    cat('Read', nrow(hsnps), 'hSNPs\n')

    cat('Joining training data..\n')
    object@gatk[hsnps, on=.(chr,pos,refnt,altnt), c('training.phgt', 'training.hap1', 'training.hap2', 'training.site') := list(i.phgt, i.hap1, i.hap2, TRUE)]
    object@gatk[is.na(training.site), training.site := FALSE]
    # index (not key) the data.table so that selecting training sites is fast
    setindex(object@gatk, training.site)
    object@training.data <- data.frame(sites=nrow(hsnps), path=path)
    object
})
