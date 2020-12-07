setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClass("SCAN2", slots=c(
    single.cell='character',
    bulk='character',
    gatk="null.or.df",
    gatk.lowmq="null.or.df",
    ab.estimates='null.or.df',
    mut.models='null.or.df',
    cigar.data='null.or.df',
    cigar.training='null.or.df',
    static.filters='null.or.df',
    static.filter='logical',
    static.filter.params='null.or.list'))

make.scan <- function(single.cell, bulk) {
    new("SCAN2", single.cell=single.cell, bulk=bulk,
        gatk=NULL, gatk.lowmq=NULL, ab.estimates=NULL, mut.models=NULL,
        cigar.data=NULL, cigar.training=NULL,
        static.filters=NULL, static.filter=NA, static.filter.params=NULL)
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

    if (!is.null(object@ab.estimates)) {
        if (!all(colnames(object@ab.estimates == c('ab', 'gp.mu', 'gp.sd'))))
            return("@ab.estimates is improperly formatted")
    }

    if (!is.null(object@mut.models)) {
        if (!all(colnames(object@mut.models == c('abc.pv', 'lysis.pv', 'lysis.beta',
            'mda.pv', 'mda.beta'))))
            return("@mut.models is improperly formatted")
    }

    # XXX: add validator for static.filters
    return(TRUE)
})


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
        cat("\n#      ", nrow(object@cigar.data), "sites\n")
        cat("#      ", nrow(object@cigar.training), "training sites\n")
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
    if (is.null(object@gatk))
        stop("must import GATK read counts first (see: read.gatk())")

    lowmq <- read.gatk.table.2sample(path, object@single.cell, object@bulk, nrows)
    lowmq <- merge(object@gatk[,c('chr', 'pos', 'refnt', 'altnt')], lowmq, all.x=TRUE)
    rownames(lowmq) <- rownames(object@gatk)

    object@gatk.lowmq <- lowmq[,c(9:10,12:13)]
    object
})


setGeneric("add.ab.estimates", function(object, path)
    standardGeneric("add.ab.estimates"))
setMethod("add.ab.estimates", "SCAN2", function(object, path) {
    if (is.null(object@gatk))
        stop("must import GATK read counts first (see: read.gatk())")

    ab <- get(load(path))
    ab <- merge(object@gatk[,c('chr','pos','refnt','altnt')], ab, all.x=TRUE)
    rownames(ab) <- rownames(object@gatk)

    # choose the AB nearest to the AF of each candidate
    # af can be NA if the site has 0 depth
    ab$gp.mu <- ifelse(!is.na(object@gatk$af) & object@gatk$af < 1/2,
        -abs(ab$gp.mu), abs(ab$gp.mu))
    ab$ab <- 1/(1+exp(-ab$gp.mu))

    object@ab.estimates <- ab[,c('ab', 'gp.mu','gp.sd')]
    object
})


setGeneric("compute.models", function(object)
    standardGeneric("compute.models"))
setMethod("compute.models", "SCAN2", function(object) {
    if (is.null(object@gatk))
        stop("must import GATK read counts first (see: read.gatk())")
    if (is.null(object@ab.estimates))
        stop("must attach allele balance estimates first (see: add.ab.estimates())")

    object@mut.models <- compute.pvs.and.betas(
        object@gatk$scalt, object@gatk$dp,
        object@ab.estimates$gp.mu, object@ab.estimates$gp.sd)
    object
})


setGeneric("add.cigar.data", function(object, sc.cigars, bulk.cigars, cigar.training)
    standardGeneric("add.cigar.data"))
setMethod("add.cigar.data", "SCAN2", function(object, sc.cigars, bulk.cigars, cigar.training) {
    if (is.null(object@gatk))
        stop("must import GATK read counts first (see: read.gatk())")

    sc <- merge(object@gatk[,c('chr','pos')], sc.cigars,
        by=c('chr', 'pos'), all.x=T)
    bulk <- merge(object@gatk[,c('chr','pos')], bulk.cigars,
        by=c('chr', 'pos'), all.x=T)
    colnames(bulk) <- paste0(colnames(bulk), '.bulk')

    # XXX: this is saved in the cigar data object for some reason
    cigar.emp.score <- function (training, test, which = c("id", "hs")) {
        xt <- training[, paste0(which, ".score.x")]
        yt <- training[, paste0(which, ".score.y")]
        x <- test[, paste0(which, ".score.x")]
        y <- test[, paste0(which, ".score.y")]
        mapply(function(xi, yi) mean(xt >= xi & yt >= yi, na.rm = T), x, y)
    }
    cigar.data <- cbind(sc[,-(1:2)], bulk[,-(1:2)])
    cigar.data$id.score.y <- cigar.data$ID.cigars / cigar.data$dp.cigars
    cigar.data$id.score.x <- cigar.data$ID.cigars.bulk / cigar.data$dp.cigars.bulk
    cigar.data$id.score <- cigar.emp.score(training=cigar.training, test=cigar.data, which='id')
    cigar.data$hs.score.y <- cigar.data$HS.cigars / cigar.data$dp.cigars
    cigar.data$hs.score.x <- cigar.data$HS.cigars.bulk / cigar.data$dp.cigars.bulk
    cigar.data$hs.score <- cigar.emp.score(training=cigar.training, test=cigar.data, which='hs')

    object@cigar.data <- cigar.data
    object@cigar.training <- cigar.training
    object
})


setGeneric("add.static.filters", function(object, min.sc.alt=2, min.sc.dp=6,
    max.bulk.alt=0, min.bulk.dp=11, exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
        standardGeneric("add.static.filters"))
setMethod("add.static.filters", "SCAN2",
function(object, min.sc.alt=2, min.sc.dp=6, max.bulk.alt=0, min.bulk.dp=11,
    exclude.dbsnp=TRUE, cg.id.q=0.05, cg.hs.q=0.05)
{
    if (is.null(object@gatk))
        stop("must import GATK read counts first (see: read.gatk())")

    object@static.filters <- data.frame(
        cigar.id.test=object@cigar.data$id.score >
            quantile(object@cigar.training$id.score, prob=cg.id.q, na.rm=T),
        cigar.hs.test=object@cigar.data$hs.score >
            quantile(object@cigar.training$hs.score, prob=cg.hs.q, na.rm=T),
        lowmq.test=is.na(object@gatk.lowmq$balt) | object@gatk.lowmq$balt <= max.bulk.alt,
        dp.test=object@gatk$dp >= min.sc.dp & object@gatk$bulk.dp >= min.bulk.dp,
        abc.test=object@mut.models$abc.pv > 0.05,
        min.sc.alt.test=object@gatk$scalt >= min.sc.alt,
        max.bulk.alt.test=object@gatk$balt <= max.bulk.alt
    )
    object@static.filter <- rowSums(object@static.filters) == 0
    object@static.filter.params <- list(
        min.sc.alt=min.sc.alt, min.sc.dp=min.sc.dp,
        max.bulk.alt=max.bulk.alt, min.bulk.dp=min.bulk.dp,
        exclude.dbsnp=exclude.dbsnp, cg.id.q=cg.id.q, cg.hs.q=cg.hs.q)
    object
})


# XXX: TODO
# make another function to add mutation type info
#    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
