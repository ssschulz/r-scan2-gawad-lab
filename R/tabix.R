# Read the first line and strip the leading '#', if it's there.
read.tabix.header <- function(tf) {
     sub('^#', '', Rsamtools::headerTabix(tf)$header)
}

# Returns a data.table
# region can only be a GRanges object with a single interval for the moment
# (we just don't have any other use cases currently).
read.tabix.data <- function(path, region, tf,
    header=read.tabix.header(tf), quiet=TRUE, ...)
{
    if ((missing(tf) & missing(path)) |(!missing(tf) & !missing(path)))
        stop('exactly one of "tf" or "path" must be specified')

    if (missing(tf)) {
        tf <- Rsamtools::TabixFile(path)
        open(tf)
    }

    if (is.null(region))
        data <- Rsamtools::scanTabix(tf)[[1]]
    else
        data <- Rsamtools::scanTabix(tf, param=region)[[1]]

    # If there are no lines in the tabix file corresponding to region (or jsut
    # no lines at all), then scanTabix returns character(0). fread() doesn't
    # like this, so replace it with an empty string.
    if (length(data) == 0)
        data <- ''

    ret <- data.table::fread(text=c(header, data), ...)

    # Only close tf if this function opened it; otherwise caller is responsible
    if (missing(tf))
        close(tf)

    if (!quiet) cat("Read", nrow(gatk), 'lines\n')

    ret
}
