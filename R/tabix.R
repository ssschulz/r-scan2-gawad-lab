# Read the first line and strip the leading '#', if it's there.
read.tabix.header <- function(tf) {
     sub('^#', '', Rsamtools::headerTabix(tf)$header)
}

# Returns a data.table
# region can only be a GRanges object with a single interval for the moment
# (we just don't have any other use cases currently).
read.tabix.data <- function(path, tf, region=NULL,
    header=read.tabix.header(tf), quiet=TRUE, ...)
{
    if ((missing(tf) & missing(path)) |(!missing(tf) & !missing(path)))
        stop('exactly one of "tf" or "path" must be specified')

    if (missing(tf)) {
        tf <- Rsamtools::TabixFile(path)
        open(tf)
        header <- read.tabix.header(tf)  # because tf didn't exist
    }

    if (is.null(region))
        data <- Rsamtools::scanTabix(tf)[[1]]
    else {
        # Important: if a chromosome is requested that isn't in the Tabix file,
        # then instead of returning empty data it throws an error. The behavior
        # we'd prefer is to return a 0-row table with the same format that would
        # otherwise be returned.
        chrs.in.file <- Rsamtools::seqnamesTabix(tf)
        chrs.in.region <- unique(seqnames(region))
        if (!all(chrs.in.region %in% chrs.in.file)) {
            # XXX: This (and many other things that interact with region) only
            # works for a single genomic interval.
            if (length(chrs.in.region) > 1)
                stop('read.tabix.data: multi-interval regions are not robustly supported. Some chromosomes in the requested region do not exist in the tabix file')
            else {
                data <- c()
            }
        } else {
            # otherwise, all regions were in the tabix file. go ahead with reading
            data <- Rsamtools::scanTabix(tf, param=region)[[1]]
        }
    }

    # If there are no lines in the tabix file corresponding to region (or jsut
    # no lines at all), then scanTabix returns character(0). fread() doesn't
    # like this, so replace it with an empty string.
    if (length(data) == 0)
        data <- ''

    ret <- data.table::fread(text=c(header, data), ...)

    # Only close tf if this function opened it; otherwise caller is responsible
    if (missing(tf))
        close(tf)

    if (!quiet) cat("Read", nrow(ret), 'lines\n')

    ret
}
