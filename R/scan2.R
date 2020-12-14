setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.logical', c('NULL', 'logical'))
setClass("SCAN2", slots=c(
    single.cell='character',
    bulk='character',
    gatk="null.or.df",
    gatk.lowmq="null.or.df",
    training.data="null.or.df",
    resampled.training.data="null.or.list",
    ab.fits='null.or.df',
    ab.estimates='null.or.df',
    mut.models='null.or.df',
    cigar.data='null.or.df',
    static.filters='null.or.df',
    static.filter='null.or.logical',
    static.filter.params='null.or.list',
    fdr.priors='null.or.df',
    fdrs='null.or.df'))


make.scan <- function(single.cell, bulk) {
    new("SCAN2", single.cell=single.cell, bulk=bulk,
        gatk=NULL, gatk.lowmq=NULL, training.data=NULL,
        ab.fits=NULL, ab.estimates=NULL, mut.models=NULL,
        cigar.data=NULL,
        static.filters=NULL, static.filter=NULL, static.filter.params=NULL,
        fdr.priors=NULL, fdrs=NULL)
}


setValidity("SCAN2", function(object) {
    if (length(object@single.cell) != 1)
        return("must provide exactly one single cell sample name")

    if (length(object@bulk) != 1)
        return("must provide exactly one bulk sample name")

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
        if (!all(colnames(object@gatk.lowmq) ==
            c('scref', 'scalt', 'bref', 'balt')))
            return("@gatk.lowmq is improperly formatted")
    }

    if (!is.null(object@training.data)) {
        if (!all(colnames(object@training.data) ==
            c('chr', 'pos', 'refnt', 'altnt', 'gt', 'hap1', 'hap2', 'dp', 'phgt')))
            return("@training.data is improperly formatted")
    }

    if (!is.null(object@ab.fits)) {
        # XXX: should also check that relevant chromosomes are present
        if (!all(colnames(object@ab.fits == c('a', 'b', 'c', 'd', 'logp'))))
            return("@ab.estimates is improperly formatted")
    }

    if (!is.null(object@ab.estimates)) {
        if (!all(colnames(object@ab.estimates == c('ab', 'gp.mu', 'gp.sd'))))
            return("@ab.estimates is improperly formatted")
    }

    if (!is.null(object@mut.models)) {
        if (!all(colnames(object@mut.models == c('abc.pv', 'lysis.pv', 'lysis.beta',
            'mda.pv', 'mda.beta'))))
            return("@mut.models is improperly formatted")
    }

    # XXX: add validator for rest of slots
    return(TRUE)
})


# various steps in the pipeline require certain slots to be filled.
# this function will ensure all slots in 'slots: (character)' are not
# NULL and will generate an error message otherwise.
check.slots <- function(object, slots) {
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
            if (s == 'ab.estimates')
                cat("must import allele balance model estimates first (see: add.ab.estimates())\n")
            if (s == 'cigar.data' | s == 'cigar.training')
                cat("must import CIGAR data first (see: add.cigar.data())\n")
            if (s == 'static.filter' | s == 'static.filters' | s == 'static.filter.params')
                cat("must apply static site filters first (see: add.static.filters())\n")
            if (s == 'fdr.priors')
                cat("must compute FDR priors first (see: compute.fdr.priors())\n")
        }
    }
    if (error.occurred)
        stop("One or more required slots are missing. See above for details.")
}

# Some getters
setGeneric("chrom", function(object) standardGeneric("chrom"))
setMethod("chrom", "SCAN2", function(object) object@gatk$chr)
setGeneric("pos", function(object) standardGeneric("pos"))
setMethod("pos", "SCAN2", function(object) object@gatk$pos)


setMethod("show", "SCAN2", function(object) {
    cat("#", is(object)[[1]], "\n")
    cat("#   Single cell ID:", object@single.cell, "\n")
    cat("#   Bulk ID:", object@bulk, "\n")
    cat("#   GATK:")
    if (is.null(object@gatk)) {
        cat(" (no data)\n")
    } else {
        cat('', nrow(object@gatk), "raw sites\n")
    }

    cat("#   GATK with low mapping quality:")
    if (is.null(object@gatk.lowmq)) {
        cat(" (no data)\n")
    } else {
        cat('', nrow(object@gatk.lowmq), "raw sites\n")
    }

    cat("#   AB model training hSNPs:")
    if (is.null(object@training.data)) {
        cat(" (no data)\n")
    } else {
        cat('', nrow(object@training.data),
            sprintf("phased sites (hap1=%d, hap2=%d)",
                sum(object@training.data$phgt=='1|0'),
                sum(object@training.data$phgt=='0|1')),"\n")
        zzz <- do.call(rbind, lapply(split(object@training.data, object@training.data$chr),
            function(td) {
                td$af <- td$hap1/td$dp
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
            cat('#        ', sum(object@resampled.training.data$selection$keep),
                'resampled hSNPs\n')
        }
    }

    cat("#   Allele balance:")
    if (is.null(object@ab.estimates)) {
        cat(" (no data)\n")
    } else {
        s <- summary(object@ab.estimates$gp.sd)
        cat('\n#       mean (0 is neutral):',
            round(mean(object@ab.estimates$gp.mu), 3), '\n')
        cat('#       uncertainty (Q25, median, Q75):',
            round(s['1st Qu.'], 3),
            round(s['Median'], 3),
            round(s['3rd Qu.'], 3), '\n')
        if (!is.null(object@training.data)) {
            cat('#       mean at training hSNPs:',
                round(mean(object@ab.estimates$gp.mu[object@gatk$training.site]), 3),
                '\n')
            cat('#       correlation with VAF at training hSNPs',
                round(cor(object@gatk$af[object@gatk$training.site],
                    object@ab.estimates$ab[object@gatk$training.site],
                    use='complete.obs'), 3), '\n')
        }
    }

    cat("#   Mutation models:")
    if (is.null(object@mut.models)) {
        cat(" (no data)\n")
    } else {
        cat(" computed\n")
    }

    cat("#   CIGAR data:")
    if (is.null(object@cigar.data)) {
        cat(" (no data)\n")
    } else {
        cat('', nrow(object@cigar.data), "sites\n")
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


internal.subset <- function(x, i) {
    if (!is.null(x@gatk))
        x@gatk <- x@gatk[i,]
    if (!is.null(x@gatk.lowmq))
        x@gatk.lowmq <- x@gatk.lowmq[i,]
    if (!is.null(x@ab.estimates))
        x@ab.estimates<- x@ab.estimates[i,]
    if (!is.null(x@mut.models))
        x@mut.models <- x@mut.models[i,]
    if (!is.null(x@cigar.data))
        x@cigar.data <- x@cigar.data[i,]
    if (!is.null(x@static.filters))
        x@static.filters <- x@static.filters[i,]
    if (!is.null(x@static.filter))
        x@static.filter <- x@static.filter[i]

    x
}


concat2 <- function(x, y) {
    if (is(x) != 'SCAN2')
        stop("x must be a SCAN2 instance")
    if (is(y) != 'SCAN2')
        stop("y must be a SCAN2 instance")
    if (x@single.cell != y@single.cell)
        stop("concat(): can only combine SCAN2 objects from the same single cell")
    if (x@bulk != y@bulk)
        stop("concat(): can only combine SCAN2 objects from the same bulk")

    x@gatk <- rbind(x@gatk, y@gatk)
    x@gatk.lowmq <- rbind(x@gatk.lowmq, y@gatk.lowmq)
    x@ab.estimates<- rbind(x@ab.estimates, y@ab.estimates)
    x@mut.models <- rbind(x@mut.models, y@mut.models)
    x@cigar.data <- rbind(x@cigar.data, y@cigar.data)
    x@static.filters <- rbind(x@static.filters, y@static.filters)
    x@static.filter <- c(x@static.filter, y@static.filter)

    # IMPORTANT: training data is never subsetted or combined
    x
}


setGeneric("concat", function(...) standardGeneric("concat"))
setMethod("concat", signature="SCAN2", function(...) {
    args <- list(...)

    if (length(args) == 0)
        stop('tried to concat() 0 objects')

    if (length(args) == 1)
        return(args[[1]])

    ret <- args[[1]]
    for (i in 2:length(args))
        ret <- concat2(ret, args[[i]])

    ret
})


# Convenience function for converting to data.frame and stitching together
# columns from all available annotations.
setGeneric("df", function(object) standardGeneric("df"))
setMethod("df", signature=c("SCAN2"), function(object) {
    check.slots(object, 'gatk')

    possible.slots <- c('gatk', 'gatk.lowmq', 'ab.estimates',
        'mut.models', 'cigar.data', 'static.filter',
        'fdr.priors', 'fdrs')

    slots <- lapply(possible.slots, function(sl) {
        # special handling to preserve name
        if (sl == 'static.filter' & !is.null(object@static.filter))
            data.frame(static.filter=object@static.filter)
        else
            slot(object=object, sl)
    })
    do.call(cbind, Filter(function(x) !is.null(x), slots))
})


setMethod("[", signature=c("SCAN2", 'ANY', 'missing', 'ANY'), function(x, i, j, drop=TRUE) {
    if (!(mode(i) %in% c('logical', 'numeric', 'integer')))
        stop("subset [ requires logical, integer indeces")
    internal.subset(x, i)
})


# Some extra work to make sure we only read in the part of the
# table relevant to these two samples. Otherwise, memory can become
# an issue for large projects.
read.gatk.table.2sample <- function(path, sc.sample, bulk.sample, nrows=-1) {
    cat("Importing GATK table..\n")

    # Step 1: read a few rows just to get column names
    gatk <- read.table(path, header=T, stringsAsFactors=F,
        colClasses=c(chr='character'), nrow=10, check.names=FALSE)
    tot.cols <- ncol(gatk)
    sc.idx <- which(colnames(gatk) == sc.sample)
    bulk.idx <- which(colnames(gatk) == bulk.sample)
    cat("Selecting columns:\n")
    for (i in 1:ncol(gatk)) {
        if (i <= 7) {
            cat(sprintf("    (%d)", i), colnames(gatk)[i], '\n')
        } else if (any(i == sc.idx + 0:2)) {
            cat(sprintf("    (%d)", i), colnames(gatk)[i], '[single cell]\n')
        } else if (any(i == bulk.idx + 0:2)) {
            cat(sprintf("    (%d)", i), colnames(gatk)[i], '[bulk]\n')
        }
    }
    
    # Step 2: really read the tables in, but only the relevant columns
    cols.to.read <- rep("NULL", tot.cols)
    # First 7 are chr, pos, dbsnp, refnt, altnt, mq, mqrs
    cols.to.read[1:7] <- c('character', 'integer', rep('character', 3), 'numeric', 'numeric')
    # Read 3 columns for the single cell, 3 columns for bulk
    cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
    cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')
    gatk <- read.table(path, header=T, stringsAsFactors=F,
        colClasses=cols.to.read, check.names=FALSE, nrows=nrows)
    cat("Read", nrow(gatk), 'lines\n')
    new.sc.idx <- which(colnames(gatk) == sc.sample)
    new.bulk.idx <- which(colnames(gatk) == bulk.sample)

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    gatk <- gatk[,c(1:7, new.sc.idx+0:2, new.bulk.idx+0:2)]
    colnames(gatk)[9:10] <- c('scref', 'scalt')
    colnames(gatk)[12:13] <- c('bref', 'balt')

    gatk
}


setGeneric("read.gatk", function(object, path,  nrows=-1)
    standardGeneric("read.gatk"))
setMethod("read.gatk", "SCAN2", function(object, path, nrows=-1) {
    gatk <- read.gatk.table.2sample(path, object@single.cell, object@bulk, nrows)

    # Add some convenient calculations
    gatk$dp <- gatk$scalt + gatk$scref
    gatk$af <- gatk$scalt / gatk$dp
    gatk$bulk.dp <- gatk$balt + gatk$bref
    gatk$bulk.af <- gatk$balt / gatk$bulk.dp
    id <- paste(gatk$chr, gatk$pos, gatk$refnt, gatk$altnt)
    rownames(gatk) <- id

    object@gatk <- gatk
    object
})

setGeneric("read.gatk.lowmq", function(object, path, nrows=-1)
    standardGeneric("read.gatk.lowmq"))
setMethod("read.gatk.lowmq", "SCAN2",
function(object, path, nrows=-1) {
    check.slots(object, 'gatk')

    lowmq <- read.gatk.table.2sample(path, object@single.cell, object@bulk, nrows)
    lowmq <- plyr::join(object@gatk[,c('chr', 'pos', 'refnt', 'altnt')], lowmq)
    rownames(lowmq) <- rownames(object@gatk)

    object@gatk.lowmq <- lowmq[,c(9:10,12:13)]
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


# Add estimates from a precomputed RDA
setGeneric("add.ab.estimates", function(object, path)
    standardGeneric("add.ab.estimates"))
setMethod("add.ab.estimates", "SCAN2", function(object, path) {
    check.slots(object, 'gatk')

    ab <- get(load(path))
    # order preservation of plyr::join is critical
    ab <- plyr::join(object@gatk[,c('chr','pos','refnt','altnt')], ab)
    rownames(ab) <- rownames(object@gatk)

    # choose the AB nearest to the AF of each candidate
    # af can be NA if the site has 0 depth
    ab$gp.mu <- ifelse(!is.na(object@gatk$af) & object@gatk$af < 1/2,
        -abs(ab$gp.mu), abs(ab$gp.mu))
    ab$ab <- 1/(1+exp(-ab$gp.mu))

    object@ab.estimates <- ab[,c('ab', 'gp.mu','gp.sd')]
    object
})


# Actually calculate the AB estimates.
# Currently, computes AB at ALL SITES.
setGeneric("compute.ab.estimates", function(object, n.cores=1)
    standardGeneric("compute.ab.estimates"))
setMethod("compute.ab.estimates", "SCAN2", function(object, n.cores=1) {
    check.slots(object, c('gatk', 'training.data', 'ab.fits'))

    sites <- object@gatk[,c('chr','pos','refnt','altnt')][1:5000,]

    # Splitting by chromosome is not for parallelization; the AB model
    # is fit separately for each chromosome and thus applies different
    # parameter estimates.
    ab <- do.call(rbind, lapply(unique(sites$chr), function(chrom) {
        hsnps <- object@training.data[object@training.data$chr == chrom,]
        fit <- object@ab.fits[chrom,,drop=FALSE]
        print(fit)

        cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
            nrow(sites), chrom, min(sites$pos), max(sites$pos)))
        time.elapsed <- system.time(z <- infer.gp1(ssnvs=sites, fit=fit,
            hsnps=hsnps, flank=1e5, verbose=TRUE,
            n.cores=n.cores))
        print(time.elapsed)
        cbind(sites, z)
    }))

    object@ab.estimates <- ab
    object
})


setGeneric("compute.models", function(object)
    standardGeneric("compute.models"))
setMethod("compute.models", "SCAN2", function(object) {
    check.slots(object, c('gatk', 'ab.estimates'))

    object@mut.models <- compute.pvs.and.betas(
        object@gatk$scalt, object@gatk$dp,
        object@ab.estimates$gp.mu, object@ab.estimates$gp.sd)
    object
})


setGeneric("compute.fdr.priors", function(object, mode='legacy')
    standardGeneric("compute.fdr.priors"))
setMethod("compute.fdr.priors", "SCAN2", function(object, mode) {
    check.slots(object, c('gatk', 'training.data', 'static.filter'))
    cat("*** FIXME: filtering out static.filter=NA; should not contain NAs.\n")
    if (mode == 'legacy') {
        # in legacy mode, all candidate sites passing a small set of pre-genotyping
        # crtieria were used.
        cand <- object@gatk[
            object@gatk$balt == 0 &
            object@gatk[,11] == '0/0' &
            object@gatk$dbsnp == '.' &
            object@gatk$scalt >= gt@static.filter.params$min.sc.alt &
            object@gatk$dp >= gt@static.filter.params$min.sc.dp &
            object@gatk$bulk.dp >= gt@static.filter.params$min.bulk.dp &
            (is.na(object@gatk.lowmq$balt) | object@gatk.lowmq$balt == 0),]
        hsnps=object@gatk[
            object@gatk$training.site &
            object@gatk$scalt >= gt@static.filter.params$min.sc.alt,]
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
            object@gatk$scalt >= gt@static.filter.params$min.sc.alt &
            object@gatk$dp >= gt@static.filter.params$min.sc.dp &
            object@gatk$bulk.dp >= gt@static.filter.params$min.bulk.dp &
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
    xt <- training[, paste0(which, ".score.x")]
    yt <- training[, paste0(which, ".score.y")]
    x <- test[, paste0(which, ".score.x")]
    y <- test[, paste0(which, ".score.y")]
    pbapply::pbmapply(function(xi, yi, bulkdp, dp)
        ifelse(dp == 0 | bulkdp == 0, 0,
            mean(xt >= xi & yt >= yi, na.rm = T)),
        x, y, test$dp.cigars.bulk, test$dp.cigars)
}


compute.cigar.scores <- function(cigar.data) {
    data.frame(
        cigar.data,
        id.score.y=cigar.data$ID.cigars / cigar.data$dp.cigars,
        id.score.x=cigar.data$ID.cigars.bulk / cigar.data$dp.cigars.bulk,
        hs.score.y=cigar.data$HS.cigars / cigar.data$dp.cigars,
        hs.score.x=cigar.data$HS.cigars.bulk / cigar.data$dp.cigars.bulk
    )
}


compute.excess.cigar <- function(cigar.data, cigar.training) {
    cat("excess indel (I/D) ops..\n")
    cigar.data$id.score <-
        cigar.emp.score(training=cigar.training, test=cigar.data, which='id')
    cat("excess clipping (H/S) ops..\n")
    cigar.data$hs.score <-
        cigar.emp.score(training=cigar.training, test=cigar.data, which='hs')
    cigar.data
}


setGeneric("add.cigar.data", function(object, sc.cigars, bulk.cigars)
    standardGeneric("add.cigar.data"))
setMethod("add.cigar.data", "SCAN2", function(object, sc.cigars, bulk.cigars) {
    check.slots(object, c('gatk', 'training.data', 'resampled.training.data'))

    # These merges just subset sc.cigars and bulk.cigars
    cat('joining CIGAR data to GATK sites..\n')
    gsc <- plyr::join(object@gatk[,c('chr','pos')], sc.cigars,
        by=c('chr', 'pos'))
    gbulk <- plyr::join(object@gatk[,c('chr','pos')], bulk.cigars,
        by=c('chr', 'pos'))
    colnames(gbulk) <- paste0(colnames(gbulk), '.bulk')

    cat('joining CIGAR data to training sites..\n')
    tsc <- plyr::join(object@training.data[,c('chr','pos')], sc.cigars,
        by=c('chr', 'pos'))
    tbulk <- plyr::join(object@training.data[,c('chr','pos')], bulk.cigars,
        by=c('chr', 'pos'))
    colnames(tbulk) <- paste0(colnames(gbulk), '.bulk')

    cat('computing CIGAR op rates for GATK sites..\n')
    gatk.cigar.data <- compute.cigar.scores(cbind(gsc[,-(1:2)], gbulk[,-(1:2)]))

    # Legacy calling used only downsampled hSNPs
    # N.B., even though every site is now scored, this computation is too
    # slow if CIGAR data from all 2-4 million sites is used.
    cat('computing CIGAR op rates for training sites..\n')
    training.cigar.data <- compute.cigar.scores(
        cbind(tsc[,-(1:2)], tbulk[,-(1:2)])[object@resampled.training.data$selection$keep,])

    cat('scoring excess CIGAR ops on all sites, using',
        nrow(training.cigar.data),
        'training sites as a reference distribution..\n')
    object@cigar.data <- compute.excess.cigar(
        cigar.data=gatk.cigar.data,
        cigar.training=training.cigar.data)
    object
})


setGeneric("add.static.filters", function(object, min.sc.alt=2, min.sc.dp=6,
    max.bulk.alt=0, min.bulk.dp=11, exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
        standardGeneric("add.static.filters"))
setMethod("add.static.filters", "SCAN2",
function(object, min.sc.alt=2, min.sc.dp=6, max.bulk.alt=0, min.bulk.dp=11,
    exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
{
    check.slots(object, c('gatk', 'cigar.data', 'gatk.lowmq', 'mut.models'))

    object@static.filters <- data.frame(
        cigar.id.test=object@cigar.data$id.score >
            quantile(object@cigar.data$id.score[object@gatk$training.site],
                prob=cg.id.q, na.rm=T),
        cigar.hs.test=object@cigar.data$hs.score >
            quantile(object@cigar.data$hs.score[object@gatk$training.site],
                prob=cg.hs.q, na.rm=T),
        lowmq.test=is.na(object@gatk.lowmq$balt) | object@gatk.lowmq$balt <= max.bulk.alt,
        dp.test=object@gatk$dp >= min.sc.dp & object@gatk$bulk.dp >= min.bulk.dp,
        abc.test=object@mut.models$abc.pv > 0.05,
        min.sc.alt.test=object@gatk$scalt >= min.sc.alt,
        max.bulk.alt.test=object@gatk$balt <= max.bulk.alt,
        dbsnp.test=!exclude.dbsnp | object@gatk$dbsnp == '.'
    )
    object@static.filter <- rowSums(object@static.filters) == ncol(object@static.filters)
    object@static.filter.params <- list(
        min.sc.alt=min.sc.alt, min.sc.dp=min.sc.dp,
        max.bulk.alt=max.bulk.alt, min.bulk.dp=min.bulk.dp,
        exclude.dbsnp=exclude.dbsnp, cg.id.q=cg.id.q, cg.hs.q=cg.hs.q)
    object
})


setGeneric("resample.training.data", function(object, M=20) 
        standardGeneric("resample.training.data"))
setMethod("resample.training.data", "SCAN2", function(object, M=20) {
    check.slots(object, c('gatk', 'training.data'))

    ret <- resample.hsnps(sites=object@gatk[object@gatk$scalt >= 2 &
            object@gatk$dp >= 6 &
            object@gatk$bulk.dp >= 11 &
            object@gatk$balt == 0,],
        hsnps=object@training.data, M=M)

    object@resampled.training.data <-
        c(ret, training.data=list(object@training.data[ret$selection$keep,]))
    object
})


setGeneric("add.training.data", function(object, path) 
        standardGeneric("add.training.data"))
setMethod("add.training.data", "SCAN2", function(object, path) {
    cat('loading', path, '\n')
    object@training.data <- get(load(path))
    cat('assigning IDs to training sites..\n')
    print(system.time(object@training.data$id <- paste(
        object@training.data$chr,
        object@training.data$pos,
        object@training.data$refnt,
        object@training.data$altnt)))

    cat('annotating GATK table\n')
    join.cols <- c('chr','pos','refnt','altnt')
    newdf <- object@training.data[,join.cols]
    newdf$training.site <- TRUE
    newdf <- plyr::join(object@gatk[,join.cols], newdf, by=join.cols)
    newdf$training.site[is.na(newdf$training.site)] <- FALSE
    object@gatk$training.site <- newdf$training.site
    object
})


# XXX: TODO
# make another function to add mutation type info
#    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
