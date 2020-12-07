setClassUnion('null.or.df', c('NULL', 'data.frame'))
setClass("SCAN2", slots=c(
    single.cell='character',
    bulk='character',
    gatk="null.or.df",
    gatk.lowmq="null.or.df",
    ab.estimates='null.or.df'))

make.scan <- function(single.cell, bulk) {
    new("SCAN2", single.cell=single.cell, bulk=bulk,
        gatk=NULL, gatk.lowmq=NULL, ab.estimates=NULL)
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
        TRUE
    }
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

    cat("#   Allele balance:\n")
    if (is.null(object@ab.estimates)) {
        cat(" (no data)\n")
    } else {
        s <- summary(object@ab.estimates$gp.sd)
        cat('#       mean (0 is neutral):',
            round(mean(object@ab.estimates$gp.mu), 3), '\n')
        cat('#       uncertainty (Q25, median, Q75):',
            round(s['1st Qu.'], 3),
            round(s['Median'], 3),
            round(s['3rd Qu.'], 3), '\n')
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
    object@ab.estimates <- ab[,c('gp.mu','gp.sd')]

    object
})

# make another function to add mutation type info
#    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
