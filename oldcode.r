# requires writing to bed and re-reading
# bgs is the result of build_annotations from the annotatr package
# i modified it to include expression quartile from BA9, GTEx
annotate.genes <- function(df, bgs, auto.delete=T) {
    require(annotatr)
    td <- paste0(tempdir())
    myd <- paste0(td, "/spmgr/")
    if (file.exists(myd) & auto.delete)
        unlink(myd, recursive=TRUE)
    dir.create(myd, recursive=TRUE)

    outf <- paste0(myd, '/muts.bed')
    options(scipen=10)
    write.table(cbind(paste0('chr', df$chr), df$pos-1, df$pos, stringsAsFactors=F),
        file=outf, sep='\t', quote=F, col.names=F, row.names=F)

    reg <- read_regions(outf, genome='hg19', format='bed')
    areg <- as.data.frame(annotate_regions(reg, bgs, ignore.strand=T))

    # annotate_regions returns a row for each matching transcript, meaning
    # each mutation can have more than one row. For now, just return row 1.
    areg$mutid <- paste(areg$seqnames, areg$start+1)
    a <- do.call(rbind, lapply(split(areg, areg$mutid), function(mut) mut[1,,drop=F]))

    a2 <- data.frame(
        chr=sub('chr', '', as.character(a$seqnames)),
        pos=a$start,
        gene=a$annot.symbol,
        gene.class=a$annot.type,
        expr.level=substr(a$annot.type, nchar(a$annot.type)-1, nchar(a$annot.type)),
        stringsAsFactors=FALSE)
    #return(a2)
    df <- merge(df, a2, by=c('chr','pos'))
    if (auto.delete)
        unlink(myd, recursive=TRUE)
    df
}







#############################################################################
# routines for joint calling
#############################################################################

joint.caller <- function(dfs, hdfs) {
    jdf <- jhelper(dfs)
    jhdf <- jhelper(hdfs)
    cutoffs <- get.joint.cutoffs(jhdf, n=length(dfs))
    jdf$pass <- apply.joint.cutoffs(jdf, cutoffs)
    dfs <- lapply(dfs, function(df) {
        # af > 1/dp is a roundabout way to enforce alt>1, because
        # at this point in the pipeline the column corresponding to
        # alt counts for this sample has been lost.
        df$jpass <- as.logical(jdf$pass * (df$af > 0) &
            df$dp.test & df$lowmq.test & df$af > 1/df$dp)
        df$jabc <- jdf$jabc
        df
    })
    dfs
}


jhelper <- function(dfs, gp.sd.penalty=0.1, gp.sd.cutoff=1) {
    n=length(dfs)
    cps <- 4 # columns per sample
    jdf <- do.call(cbind, lapply(dfs,function(df) df[,c('af','dp','abc.pv','gp.sd')]))

    jdf$ncells <- rowSums(!apply(jdf[,seq(1,n*cps,cps)], 1:2, is.nan) & jdf[,seq(1,n*cps,cps)]>0)

    # column reordering makes it easy to recover af, dp, pv
    jdf$jabc <- apply(jdf[,c(seq(1,n*cps,cps), seq(2,n*cps,cps), seq(3,n*cps,cps), seq(4,n*cps,cps))], 1,
        function(row) {
            afs <- row[1:n]
            dps <- row[(n+1):(2*n)]
            pvs <- row[(2*n+1):(3*n)]
            sds <- row[(3*n+1):(4*n)]
            n.penalties <- sum(sds[dps>0 & afs>0] >= gp.sd.cutoff)
            prod(pvs[dps>0 & afs>0]) * prod(rep(gp.sd.penalty, n.penalties))
    })
    jdf
}

# n - number of samples
get.joint.cutoffs <- function(jdf, n, q=0.1) {
    counts <- sapply(2:n, function(i) sum(jdf$ncells == i))
    if (any(counts) == 0) {
        cat("counts for ncells:\n")
        for (i in 2:n)
            cat(sprintf("ncells=%2d, count=%d\n", i, counts[i-1]))
        stop("Found ncell count=0 (see above).  Try increasing the number of hSNP spikeins")
    }

    cutoffs <- sapply(2:n, function(i)
        quantile(jdf$jabc[jdf$ncells == i], prob=q))
    # ncells=0,1: don't pass any of these. they will be determined only by
    # single sample statistics
    data.frame(ncells=0:n, cutoff=c(Inf, Inf, cutoffs))
}


apply.joint.cutoffs <- function(jdf, cutoffs) {
    mapply(function(ncells, jabc)
        jabc >= cutoffs$cutoff[cutoffs$ncells == ncells],
        ncells=jdf$ncells, jabc=jdf$jabc)
}
