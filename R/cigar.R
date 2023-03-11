# Functions for dealing with the CIGAR I/D and H/S filters.
#
# In general, these filters are designed to test whether the single cell
# has an excessive number of indel (I/D) or clipping (H/S) CIGAR operations
# compared to the matched bulk.
#
# N.B. mismatches, as would be caused by a somatic SNV, generate an M CIGAR
# op, so the bulk and single cell should not have different CIGAR profiles.
# However, a somatic indel *does* generate an I or D OP that would not be
# seen in the bulk. So this filter is not well-designed for indels; yet it
# is still applied.
#
# Currently, several CIGAR scoring operations cannot be applied to chunked
# objects, so we are stuck using relatively slow calculations and small
# training sets.


cigar.emp.score <- function (training, test, which = c("id", "hs"), legacy=FALSE, quiet=FALSE) {
    xt <- training[[paste0(which, ".score.x")]]
    yt <- training[[paste0(which, ".score.y")]]
    x <- test[[paste0(which, ".score.x")]]
    y <- test[[paste0(which, ".score.y")]]

    dp <- test$dp.cigars
    bulkdp <- test$dp.cigars.bulk

    progressr::with_progress({
        if (!quiet) p <- progressr::progressor(along=1:(length(x)/100))
        # in legacy mode, NA and NaN values in x or y weren't caught
        if (legacy) {
            ret <- future.apply::future_mapply(function(dp, bulkdp, x, y, i) {
                if (!quiet & i %% 100 == 1) p()
                mean(xt >= x & yt >= y, na.rm = T)
            }, test$dp.cigars, test$dp.cigars.bulk, x, y, 1:length(x))
        } else {
            ret <- future.apply::future_mapply(function(dp, bulkdp, x, y, i) {
                if (!quiet & i %% 100 == 1) p()
                ifelse(dp == 0 | bulkdp == 0, 0,
                    mean(xt >= x & yt >= y, na.rm = T))
            }, test$dp.cigars, test$dp.cigars.bulk, x, y, 1:length(x))
        }
    }, enable=!quiet)

    # future_mapply returns list() when x is length 0. make this a 0-length
    # numeric so quantile() doesn't fail on it later.
    if (length(ret) == 0)
        return(numeric(0))
    else
        return(ret)
}


compute.cigar.scores <- function(cigar.data) {
    cigar.data[, c('id.score.y', 'id.score.x', 'hs.score.y', 'hs.score.x') :=
        list(ID.cigars / dp.cigars,
             ID.cigars.bulk / dp.cigars.bulk,
             HS.cigars / dp.cigars,
             HS.cigars.bulk / dp.cigars.bulk)]
}


read.cigar.data <- function(path, region, quiet=FALSE) {
    if (!quiet) cat('Importing CIGAR stats from', path, '\n')
    col.classes <- c('character', rep('integer', 6))
    read.tabix.data(path=path, region=region, quiet=quiet, colClasses=col.classes)
}

        
setGeneric("add.cigar.data", function(object, sc.cigars.path, bulk.cigars.path, quiet=FALSE)
    standardGeneric("add.cigar.data"))
setMethod("add.cigar.data", "SCAN2", function(object, sc.cigars.path, bulk.cigars.path, quiet=FALSE) {
    check.slots(object, 'gatk')

    sc <- read.cigar.data(sc.cigars.path, region=object@region, quiet=quiet)
    bulk <- read.cigar.data(bulk.cigars.path, region=object@region, quiet=quiet)

    if (!quiet) cat('joining CIGAR data..\n')
    object@gatk[sc, on=c('chr', 'pos'),
        c('M.cigars', 'ID.cigars', 'HS.cigars', 'other.cigars', 'dp.cigars') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]
    object@gatk[bulk, on=c('chr', 'pos'),
        c('M.cigars.bulk', 'ID.cigars.bulk', 'HS.cigars.bulk', 'other.cigars.bulk', 'dp.cigars.bulk') :=
            list(i.M.cigars, i.ID.cigars, i.HS.cigars, i.other.cigars, i.dp.cigars)]

    if (!quiet) cat('computing CIGAR op rates..\n')
    compute.cigar.scores(object@gatk)  # modifies by reference
    object@cigar.data <- data.frame(sc.sites=nrow(sc), bulk.sites=nrow(bulk),
        sc.path=sc.cigars.path, bulk.path=bulk.cigars.path)
    object
})


# It would be nice to use all training sites for this; however, computing excess
# cigar scores is O(n^2), n=number of null sites.
cigar.get.null.sites <- function(object, path=NULL, legacy=TRUE, quiet=FALSE) {
    if (is.null(path)) {
        check.slots(object, c('gatk', 'cigar.data'))
        if (legacy) {
            null.sites <- object@gatk[resampled.training.site==TRUE]
        } else {
            null.sites <- object@gatk[training.site==TRUE]
        }
    } else {
        null.sites <- data.table::fread(path)
    }

    null.sites
}


setGeneric("compute.excess.cigar.scores", function(object, path=NULL, legacy=TRUE, quiet=FALSE)
    standardGeneric("compute.excess.cigar.scores"))
setMethod("compute.excess.cigar.scores", "SCAN2", function(object, path=NULL, legacy=TRUE, quiet=FALSE) {
    check.slots(object, c('gatk', 'static.filter.params', 'cigar.data'))
    null.sites <- cigar.get.null.sites(object, path, legacy, quiet)
    compute.cigar.scores(null.sites)
    if (!quiet) {
        if (legacy) {
            cat(sprintf('LEGACY: computing CIGAR op rates only at resampled training sites (n=%d)..\n',
                nrow(null.sites)))
        } else {
            cat(sprintf('WARNING: using the full set of training het germline sites (n=%d) for CIGAR op calculations may be prohibitively slow\n',
                nrow(null.sites)))
        }
    }
    muttypes <- c('snv', 'indel')
    object@excess.cigar.scores <- setNames(lapply(muttypes, function(mt) {
        null.sites.mt <- null.sites[muttype == mt]
        pc <- perfcheck("excess CIGAR ops",
            object@gatk[muttype == mt, c('id.score', 'hs.score') := list(
                    cigar.emp.score(training=null.sites.mt, test=.SD, which='id', quiet=quiet, legacy=legacy),
                    cigar.emp.score(training=null.sites.mt, test=.SD, which='hs', quiet=quiet, legacy=legacy)
            )],
            report.mem=FALSE)
        if (!quiet) cat(pc, '\n')

        # it is only necessary to score training sites to get the
        # quantiles for test cutoff.
        null.sites.mt[, c('id.score', 'hs.score') := list(
            cigar.emp.score(training=null.sites.mt, test=null.sites.mt, which='id', quiet=quiet, legacy=legacy),
            cigar.emp.score(training=null.sites.mt, test=null.sites.mt, which='hs', quiet=quiet, legacy=legacy)
        )]
        sfp <- object@static.filter.params[[mt]]
        data.frame(sites=nrow(object@gatk[muttype==mt]),
            null.sites=nrow(null.sites.mt),
            # These must not be calculated on chunked objects; the full
            # set of null sites must be available.
            id.score.q=quantile(null.sites.mt$id.score, prob=sfp$cg.id.q, na.rm=TRUE),
            hs.score.q=quantile(null.sites.mt$hs.score, prob=sfp$cg.hs.q, na.rm=TRUE),
            legacy=legacy)
    }), muttypes)
    object
})
