# Some extra work to make sure we only read in the part of the
# table relevant to one sample. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
#
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
#
# keep.coords - whether to keep the chrom/pos columns. these files are
#    basepair resolution, so doubling the memory required for reading
#    can have a large impact on RAM needed.
read.depth.2sample <- function(path, sc.sample, bulk.sample, keep.coords=TRUE, region=NULL, quiet=FALSE) {
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
    close(tf)
    col.strings <- strsplit(header, '\t')[[1]]

    if (!(sc.sample %in% col.strings))
        stop(paste('single cell sample', sc.sample, 'not found in table'))
    sc.sample.idx <- which(col.strings == sc.sample)

    if (!(bulk.sample %in% col.strings))
        stop(paste('bulk sample', bulk.sample, 'not found in table'))
    bulk.sample.idx <- which(col.strings == bulk.sample)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (i == sc.sample.idx) {
                cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
            } else if (i == bulk.sample.idx) {
                cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
            } else {
                cat(sprintf("    (%d)\n", i))
            }
        }
    }
    
    cols.to.read <- rep("NULL", length(col.strings))
    cols.to.read[c(sc.sample.idx, bulk.sample.idx)] <- c('integer', 'integer')

    # Do we need chromosome/position?
    if (keep.coords)
        cols.to.read[1:2] <- c('character', 'integer')

    # Reading in a somewhat preparsed GATK DepthOfCoverage table
    gatk.doc <- read.tabix.data(path=path, region=region, header=header, quiet=quiet, colClasses=cols.to.read)

    colorder <- c(sc.sample, bulk.sample)
    if (keep.coords)
        colorder <- c('chr', 'pos', colorder)
    setcolorder(gatk.doc, colorder)

    gatk.doc
}


# Some extra work to make sure we only read in the part of the
# table relevant to one sample. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
#
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
digest.depth.2sample <- function(path, sc.sample, bulk.sample, clamp.dp=500, region=NULL, quiet=FALSE) {
    gatk.doc <- read.depth.2sample(path=path,
        sc.sample=sc.sample, bulk.sample=bulk.sample,
        keep.coords=FALSE, region=region, quiet=quiet)

    if (nrow(gatk.doc) > 0) {
        # Set maximum depth to clamp.dp
        gatk.doc <- gatk.doc[, lapply(.SD, pmin, clamp.dp)]

        # add points (0,0) ... (clamp.dp,clamp.dp) to the gatk depthofcoverage output
        # so that the result of R's table() will at least be (clamp.dp x clamp.dp)
        # in dimension.
        # subtracting one from the diagonal easily removes this afterward.
        gatk.doc <- rbind(gatk.doc, data.table(0:clamp.dp, 0:clamp.dp), use.names=FALSE)
    } else {
        gatk.doc <- data.table(0:clamp.dp, 0:clamp.dp)
    }

    dptab <- table(gatk.doc) - diag(clamp.dp+1)
    dptab
}


# Returns a GRanges object of all ranges passing the minimum depth cutoffs
# given in min.sc.dp and min.bulk.dp.
compute.callable.region <- function(path, sc.sample, bulk.sample, min.sc.dp, min.bulk.dp, region=NULL, quiet=FALSE)
{
    gatk.doc <- read.depth.2sample(path=path,
        sc.sample=sc.sample, bulk.sample=bulk.sample,
        keep.coords=TRUE, region=region, quiet=quiet)

    colnames(gatk.doc)[3:4] <- c('sc.dp', 'bulk.dp')
    gatk.doc <- gatk.doc[sc.dp >= min.sc.dp & bulk.dp >= min.bulk.dp]
    # GPos() doesn't support reduce() so sadly it seems we don't have the
    # option to use it for more memory efficiency.
    g <- reduce(GRanges(seqnames=gatk.doc$chr, ranges=IRanges(start=gatk.doc$pos, end=gatk.doc$pos)))
    g
}



setGeneric("mean.coverage", function(object) standardGeneric("mean.coverage"))
setMethod("mean.coverage", "SCAN2", function(object) {
    check.slots(object, 'depth.profile')

    # rows correspond to single cell depth values, starting at 0 (row 1 = #bases at dp=0,
    # row 2 = #bases at dp=1, ...). rownames() actually records dp values as strings.
    dptab <- object@depth.profile$dptab
    c(mean.coverage=sum(rowSums(dptab) * as.numeric(rownames(dptab))) / sum(dptab))
})
