# Read the header and strip the leading '#', if it's there.
# If the header is multilined (e.g., VCFs), then we assume the last line
# contains column names.
read.tabix.header <- function(tf) {
     sub('^#', '', tail(Rsamtools::headerTabix(tf)$header, 1))
}


# Rsamtools::scanTabix does not provide a mechanism to read a subset of
# columns. This causes major memory over-usage in large projects with
# 10s-100s of cells when reading joint/integrated tables.
#
# If colClasses was specified, use a system call to the tabix binary
# (meaning tabix must be installed on the system; the Rsamtools package
# likely is not sufficient) piped through cut to get only the relevant
# columns.
tabix.read.only.cols <- function(path, header, colClasses, region=NULL) {
print('got header ----')
print(header)

    if (!is.null(region) & length(region) > 1)
        stop('region must contain only a single interval')

    command <- paste('tabix', path)
    if (!is.null(region))
        command <- paste(command,
            sprintf('%s:%d-%d', seqnames(region)[1], start(region)[1], end(region)[1]))

    # handle two types of colClasses:
    #  1. character vector of types, 1 per column. If any column is not specified
    #     then it is not read.
    #  2. list of column types with names in the form of, e.g., character=c('column 1', 'column 7', ...)
    #     IMPORTANTLY: it doesn't matter what is specified in the list form; we DO NOT ALLOW
    #     NULL SKIPPING in the list form. So all columns will be read.
    if (!missing(colClasses) & !is.list(colClasses)) {
        if (is.character(colClasses)) {
            new.col.classes <- rep('NULL', length(header))
            new.col.classes[1:length(colClasses)] <- colClasses
        } else {
            stop('colClasses must either be a character vector or a list')
        }
        cut.command <- paste0('cut -f', paste(which(new.col.classes != 'NULL'), collapse=','))
        command <- paste(command, "|", cut.command)
    } 

print(command)
    system(command, intern=TRUE)
}


# Returns a data.table
# region can only be a GRanges object with a single interval for the moment
# (we just don't have any other use cases currently).
read.tabix.data <- function(path, header, region=NULL, quiet=TRUE, ...)
{
    if (missing(header)) {
        tf <- Rsamtools::TabixFile(path)
        open(tf)
        header <- read.tabix.header(tf)
        close(tf)
    }

    if (is.null(region)) {
        #data <- Rsamtools::scanTabix(tf)[[1]]
        data <- tabix.read.only.cols(path=path, header=header, region=NULL, ...)
    } else {
        # Important: if a chromosome is requested that isn't in the Tabix file,
        # then instead of returning empty data it throws an error. The behavior
        # we'd prefer is to return a 0-row table with the same format that would
        # otherwise be returned.
        tf <- Rsamtools::TabixFile(path)
        open(tf)
        chrs.in.file <- Rsamtools::seqnamesTabix(tf)
        close(tf)
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
            #data <- Rsamtools::scanTabix(tf, param=region)[[1]]
            data <- tabix.read.only.cols(path=path, header=header, region=region, ...)
        }
    }

    # If there are no lines in the tabix file corresponding to region (or jsut
    # no lines at all), then scanTabix returns character(0). fread() doesn't
    # like this, so replace it with an empty string.
    if (length(data) == 0)
        data <- ''

cat("header = ------------------------\n")
print(header)

cat("data = ------------------------\n")
str(data)

    ret <- data.table::fread(text=c(header, data), ...)

    if (!quiet) cat("Read", nrow(ret), 'lines\n')

str(ret)
    ret
}
