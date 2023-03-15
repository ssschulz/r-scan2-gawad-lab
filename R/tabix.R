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
    if (!is.null(region) & length(region) > 1)
        stop('region must contain only a single interval')

    # As far as I know, there is no way or the tabix command to return
    # the entire file. So when region is null, we can't use the same
    # tabix system call method. This is unfortunate since it'd be nice to
    # not have special behavior when region=NULL, but I suppose it can't
    # be avoided.
    # tail - don't print the header line, ASSUMING it's only 1 line
    if (is.null(region)) {
        command <- paste('gunzip -c', path, "| grep -v '^#'")
    } else {
        command <- paste('tabix', path)
        if (!is.null(region))
            command <- paste(command,
                sprintf('%s:%d-%d', seqnames(region)[1], start(region)[1], end(region)[1]))
    }

    # handle two types of colClasses:
    #  1. character vector of types, 1 per column. If any column is not specified
    #     then it is not read.
    #  2. list of column types with names in the form of, e.g., character=c('column 1', 'column 7', ...)
    #     IMPORTANTLY: it doesn't matter what is specified in the list form; we DO NOT ALLOW
    #     NULL SKIPPING in the list form. So all columns will be read.
    if (!is.null(colClasses) & !is.list(colClasses)) {
        if (is.character(colClasses)) {
            new.col.classes <- rep('NULL', length(header))
            new.col.classes[1:length(colClasses)] <- colClasses
        } else {
            stop('colClasses must either be a character vector or a list')
        }
        cut.command <- paste0('cut -f', paste(which(new.col.classes != 'NULL'), collapse=','))
        command <- paste(command, "|", cut.command)
    } 

    system(command, intern=TRUE)
}


# Returns a data.table
# region can only be a GRanges object with a single interval for the moment
# (we just don't have any other use cases currently).
read.tabix.data <- function(path, header, region=NULL, colClasses=NULL, quiet=TRUE, ...)
{
    if (missing(header)) {
        tf <- Rsamtools::TabixFile(path)
        open(tf)
        header <- read.tabix.header(tf)
        close(tf)
    }

    if (is.null(region)) {
        data <- tabix.read.only.cols(path=path, header=header, region=NULL, colClasses=colClasses) 
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
            data <- tabix.read.only.cols(path=path, header=header, region=region, colClasses=colClasses) 
        }
    }

    if (!quiet) cat("Read", length(data), 'lines\n')

    # If there are no lines in the tabix file corresponding to region (or just
    # no lines at all), then scanTabix returns character(0). fread() doesn't
    # like this, so replace it with an empty string.
    if (length(data) == 0)
        data <- ''

    # tabix.read.only.cols has already dropped any NULL colClasses
    # this function does not fully implement colClasses features as in
    # read.table and fread
    if (!is.null(colClasses)) {
        if (!is.list(colClasses)) { # Remember: our list colClasses version does not support NULL skipping
            # support the case where user only specifies a subset of colClasses
            header <- paste(strsplit(header, '\t')[[1]][1:length(colClasses)][colClasses != 'NULL'], collapse='\t')
            colClasses <- colClasses[colClasses != 'NULL']
        }
        # header=TRUE is critical: if any column names are determined by data.table to not
        # be of type character, then the first row is NOT considered a header row. E.g.,
        # if the header contains a column name that is a number (like sample ID=1234), then
        # fread will reject the entire header line and treat it as a data line.
        ret <- data.table::fread(text=c(header, data), colClasses=colClasses, header=TRUE, ...)
    } else {
        ret <- data.table::fread(text=c(header, data), header=TRUE, ...)
    }

    # There is a strange special case in data.table::fread where if there is only
    # a single column and no data (data=''), then instead of resolving to an empty,
    # one column data.table, it produces a 1 row data.table with NA as the column.
    # Handle this case by removing the row.
    if (length(data) == 1 & length(strsplit(header, '\t')[[1]]) == 1) {
        if (data == '') {
            ret <- ret[-1]
        }
    }

    # Handle a corner case that causes a site to show up in non-overlapping regions.
    # If an indel spans a region boundary, then tabix returns the indel for each region
    # spanned. E.g., an indel such as
    #     5   119999998   .   GAAAAC  G
    # starts at 119999998 and ends at 120000003. If a chunk boundary exists at 120000001,
    # then this site will be returned twice: for the chunk before and after the boundary:
    #     $ tabix mmq60.tab.gz 5:110000001-120000000|tail -n 1|cut -f1-5
    #     5   119999998   .   GAAAAC  G
    #     $ tabix mmq60.tab.gz 5:120000001-130000000|head -1|cut -f1-5
    #     5   119999998   .   GAAAAC  G
    # If we allow this, it will cause duplicates. It doesn't matter which chunk the site
    # is assigned to so long as it's consistent.
    #
    # For now, print out a notification
    if (!is.null(region)) {
        out.of.chunk.site <- data$pos < start(region) | data$pos > end(region)
        s <- sum(out.of.chunk.site)
        if (s > 0) {
            cat(paste0("*** INFO: removing the following ", s, " sites, (out of requested chunk bounds: ", seqnames(region)[1], ":", start(region)[1], "-", end(region[1]), "). this is normal behavior for tabix, but is not desirable."))
            print(data[out.of.chunk.site==TRUE])
            data <- data[out.of.chunk.site == FALSE]
        }
    }
    
    ret
}
