# Some extra work to make sure we only read in the part of the
# table relevant to one sample. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
#
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
read.depth.1sample <- function(path, sc.sample, bulk.sample, clamp.dp=500, region=NULL, quiet=FALSE) {
    # Step 1: just get the header and detect the columns corresponding to sample.id
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
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
    
    # Step 2: really read the tables in, but only the relevant columns
    cols.to.read <- rep("NULL", length(col.strings))
    # We don't care about chromosome/position here, only depth counts
    cols.to.read[c(sc.sample.idx, bulk.sample.idx)] <- c('integer', 'integer')

    # Reading in a somewhat preparsed GATK DepthOfCoverage table
    gatk.doc <- read.tabix.data(tf=tf, region=region, header=header, quiet=quiet, colClasses=cols.to.read)
    close(tf)

    # Standardize on 1st column: single cell, 2nd column: bulk
    setcolorder(gatk.doc, c(sc.sample, bulk.sample))

    # Set maximum depth to clamp.dp
    gatk.doc <- apply(gatk.doc, 2, pmin, clamp.dp)

    # add points (0,0) ... (clamp.dp,clamp.dp) to the gatk depthofcoverage output
    # so that the result of R's table() will at least be (clamp.dp x clamp.dp)
    # in dimension.
    # subtracting one from the diagonal easily removes this afterward.
    gatk.doc <- rbind(gatk.doc, data.table(0:clamp.dp, 0:clamp.dp), use.names=FALSE)

    dptab <- table(gatk.doc) - diag(clamp.dp+1)

    dptab
}
