# for 96 dimensional mut sigs
mutsig.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey', 'chartreuse3', 'pink2'), each=16)

# Returns the plotted spectrum invisibly
plot.3mer <- function(x, no.legend=FALSE, ...) {
    bases <- c("A", 'C', 'G', 'T')

    # need to make a table of all possible contexts because they may not
    # be observed after filtering.
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each=4),
                      rep(c('C', 'T'), each=48),
                      rep(bases, times=4),
                      ":",
                      rep(c('C', 'T'), each=48),
                      ">",
                      c(rep(c('A', 'G', 'T'), each=16),
                        rep(c('A', 'C', 'G'), each=16)))
    t2 <- table(x$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[,2])]
    p <- barplot(t, las=3, col=mutsig.cols, names.arg=tn[order(tn[,2]), 1], space=0.5, border=NA, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    if (!no.legend)
        legend('topright', ncol=2, legend=sort(unique(tn[,2])),
            fill=mutsig.cols[seq(1, length(mutsig.cols), 16)])
    invisible(t)
}


# Returns the ordered channel values invisibly
plot.indel <- function(iclass, tsb=F, proc, reduce.to.id83=FALSE, xaxt='n',
    col, border, make.plot=TRUE, ...) {
    # ID83 plotting order
    iclass.order <- paste(
        c(rep(1,24), rep(rep(2:5,each=6),2),c(2,3,3,4,4,4,5,5,5,5,5)),
        c(rep(c('Del', 'Ins'), each=12), rep(c('Del', 'Ins'), each=24), rep('Del', 11)),
        c(rep(rep(c('C','T'), each=6), 2), rep('R',48), rep('M', 11)),
        c(rep(0:5, 12), c(1, 1,2, 1:3, 1:5)),
        sep=':')
    # ID415 (includes transcribed strand)
    if (tsb)
        iclass.order <- as.vector(rbind(
            paste0("T:", iclass.order),
            paste0("U:", iclass.order)))

    if (reduce.to.id83)
        iclass <- substr(iclass,3,11)

    if (missing(proc)) {
        proc <- sapply(iclass.order, function(ico) sum(iclass==ico))
    } else
        proc <- proc[iclass.order]
    if (missing(border))
        border <- rep(c('#FBBD75', '#FC7F24', '#B0DB8E', '#3B9F36', '#FBC9B6', '#F8896D', '#EE453A', '#B91C22', '#C5D4E4', '#8DBAD2', '#4D98C6', '#1D65A8', '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'), c(rep(6,12), 1,2,3,5))
    if (missing(col))
        col <- rep(c('#FBBD75', '#FC7F24', '#B0DB8E', '#3B9F36', '#FBC9B6', '#F8896D', '#EE453A', '#B91C22', '#C5D4E4', '#8DBAD2', '#4D98C6', '#1D65A8', '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'), c(rep(6,12), 1,2,3,5))
    if (tsb)
        mutsig.cols <- as.vector(rbind(mutsig.cols, "#D6C2C2"))

    if (make.plot) {
        par(mar=c(4, 4, 1, 1))
        p <- barplot(proc, las = 3, col = col, names.arg = iclass.order,
            space = 0.5, border = border, cex.names=0.7, ...) #xaxt=xaxt, ...)
        abline(v = (p[seq(6, length(p) - 11, 6)] + p[seq(7, length(p)-10,6)])/2,
        col="grey")
        mtext(text=c('del C', 'del T', 'ins C', 'ins T', 'del 2', 'del 3', 'del 4',
            'del 5+', 'ins 2', 'ins 3', 'ins 4', 'ins 5+', 'microhom.'),
            side=1, at=c(mean(p[1:6]), mean(p[7:12]), mean(p[13:18]), mean(p[19:24]), mean(p[25:30]), mean(p[31:36]), mean(p[37:42]), mean(p[43:48]), mean(p[49:54]), mean(p[55:60]), mean(p[61:66]), mean(p[67:72]), mean(p[73:83])))
    }
    invisible(proc)
}

#######################################################################
# Plots for AB model and SCAN-SNV internal estimation procedures
#######################################################################

# Plot AB model and error distributions for a single site
plot.ab <- function(ab) {
    layout(matrix(1:4,nrow=2,byrow=T))
    td <- ab$td[order(ab$td$dp),]
    plot(td$dp, td$mut, type='l', ylim=range(td[,c('mut', 'err1', 'err2')]),
        xlab="Depth", ylab="Model probabiblity")
    lines(td$dp, td$err, lty='dotted', col=2)
    plot(ab$alphas, ab$betas, xlab='FP rate', ylab='Power')
    abline(h=0, lty='dotted')
    plot(ab$alphas, ab$betas, log='x', xlab='log(FP rate)', ylab='Power')
    abline(h=0, lty='dotted')
}

# Plot histograms of the numbers of estimated TPs and FPs in
# a particular subset of SNV candidates (typically stratified
# by VAF and DP).
plot.fcontrol <- function(fc) {
    layout(matrix(1:(1+length(fc$pops)), nrow=1))
    plot(fc$binmids, fc$g/sum(fc$g),
        ylim=range(c(fc$g/sum(fc$g), fc$s/sum(fc$s))),
        type='l', lwd=2)
    lines(fc$binmids, fc$s/sum(fc$s), col=2, lwd=2)

    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        barplot(names.arg=fc$binmids, t(pop), col=1:2,
            main=sprintf('Assumption: ~%d true sSNVs', sum(round(pop[,1],0))), las=2)
        legend('topright', fill=1:2, legend=c('Ntrue', 'Nartifact'))
    }
}

# using the (alpha, beta) relationships and fcontrol population
# estimations, determine average sensitivity per AF with a
# (theoretically) controlled FDR
plot.fdr <- function(fc, dps=c(10,20,30,60,100,200), target.fdr=0.1, div=2) {
    afs <- fc$binmids
    layout(matrix(1:(3*length(fc$pops)), nrow=3))
    for (i in 1:length(fc$pops)) {
        # from fcontrol: pops[[i]] has rows corresponding to afs
        pop <- fc$pops[[i]]
        l <- lapply(dps, function(dp) {
            sapply(1:length(afs), function(i)
                match.fdr(afs[i], dp, nt=pop[i,1], na=pop[i,2],
                    target.fdr=target.fdr, div=div)
            )
        })
        matplot(x=afs, sapply(l, function(ll) ll[3,]), type='l', lty=1,
            main=sprintf("Assuming %d true sSNVs", sum(pop[,1])),
            xlab="AF (binned)", ylab="FDR", ylim=c(0, 1.1*target.fdr))
        abline(h=target.fdr, lty='dotted')
        matplot(x=afs, sapply(l, function(ll) ll[1,]), type='l', lty=1,
            xlab="AF (binned)", ylab="log(alpha)", log='y', ylim=c(1e-5,1))
        abline(h=10^-(1:5), lty=2)
        matplot(x=afs, sapply(l, function(ll) ll[2,]), type='l', lty=1,
            xlab="AF (binned)", ylab="Power", ylim=0:1)
    }
}


# Plot the AB model in a genomic window with confidence bands.
# Optionally: plot a candidate/putative sSNV.
# XXX: TODO: restore the ability to compute GP mu/sd at many additional
# points in the region to give smooth bands/lines.
setGeneric("plot.region", function(object, site=NA, chrom=NA, pos=NA, upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=FALSE)
    standardGeneric("plot.region"))
setMethod("plot.region", "SCAN2", function(object, site=NA, chrom=NA, pos=NA,
    upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=FALSE)
{
    check.slots(object, c('gatk', 'training.data', 'ab.estimates'))
    if (recompute)
        check.slots(object, 'ab.fits')

    if (!missing(site) & (!missing(chrom) | !missing(pos)))
        stop("either site or (chrom,pos) must be specified, but not both")
    if ((!missing(chrom) & missing(pos)) | (missing(chrom) & !missing(pos)))
        stop("both chrom and pos must be specified")

    if (!missing(site)) {
        chrom <- object@gatk[site,]$chr
        pos <- object@gatk[site,]$pos
    }

    # Sites at which AB was estimated in the genotyper.
    # Each site in 'd' will be plotted with a point.
    # Don't keep sites with no reads in this sample, unless it's a training
    # site. The majority of these 0 read, non-germline sites are here
    # because they had reads in a different sample.
    # Also, homozygous sites can be confusing since they appear at 0 or 1.
    d <- object@gatk[chrom(object) == chrom &
        pos(object) >= pos - upstream & pos(object) <= pos + upstream &
        (object@gatk$training.site | object@gatk$scalt > 0) &
        object@gatk[,11] != '1/1',]  # column 11 is bulk GT

    # Old code that computed GP on a fine grid surrounding the target
    if (recompute) {
        cat("estimating AB in region..\n")
        # ensure that we estimate at exactly pos and at all sites in 'd'
        # other loci are not sites reported in the GATK table, they are
        # only there to make smooth lines.
        est.at <- c(seq(pos - upstream, pos-1, length.out=n.gp.points/2), pos,
                    seq(pos+1, pos + downstream, length.out=n.gp.points/2),
                    d$pos)
        est.at <- sort(unique(est.at))
        fit.chr <- object@ab.fits[chrom,]
        gp <- as.data.frame(infer.gp1(ssnvs=data.frame(pos=est.at),
            fit=fit.chr,
            hsnps=object@training.data[object@training.data$chr == chrom,],
            flank=gp.extend, max.hsnp=150))
        gp$chr <- chrom
        gp$pos <- est.at
        gp$ab <- 1/(1+exp(-gp$gp.mu))

        # Need to reflect sites to match the GP
        d <- merge(d, gp[,c('pos','ab')], all.x=TRUE)  # attach gp estimtes
        d$af <- ifelse(abs(d$af - d$ab) <= abs(1-d$af - d$ab), d$af, 1-d$af)
    } else {
        # Rely on precomputed GP mu/sd (relevant for ALLSITES mode)
        gp <- cbind(
            object@gatk[chrom(object) == chrom &
                pos(object) >= pos - upstream & pos(object) <= pos + upstream,],
            object@ab.estimates[chrom(object) == chrom &
                pos(object) >= pos - upstream & pos(object) <= pos + upstream,]
        )
    }

    plot.gp.confidence(df=gp, add=FALSE)
    points(d$pos, d$af, pch=20,
        cex=ifelse(d$training.site, 1, 1.5),
        col=2 - d$training.site, ylim=c(-0.2,1))
    # 5*max : restrict the depth to the bottom 20% of plot
    lines(d$pos, d$dp/(5*max(d$dp)), type='h', lwd=2)
    text(d$pos[which.max(d$dp)], 1/5, max(d$dp))
    abline(h=30/(5*max(d$dp)), lty='dotted', col='grey')

    lines(gp$pos, gp$ab/2, lwd=2, col=2)
    lines(gp$pos, (1-gp$ab)/2, lwd=2, col=2)

    # If a site was given, plot it specially
    if (!missing(site)) {
        abline(v=pos, lty='dotted')
        abline(h=d$af[d$pos==pos], lty='dotted')
        points(pos, d$af[d$pos == pos], pch=4, cex=1.5, lwd=2, col=3)
    }

    legend('topright', legend=c('Training hSNP', 'Other site', 'Target site'),
        pch=c(20,20,4), col=1:3, pt.cex=c(1,1.5,1.5))
})

# Add 95% confidence bands for AB model to plot.
# NOTE: the 95% probability interval is in the NORMAL space, not the
# fraction space [0,1]. After the logistic transform, the region in
# fraction space may not be an HDR.
plot.gp.confidence <- function(pos, gp.mu, gp.sd, df, sd.mult=2,
    logspace=FALSE, tube.col=rgb(0.9, 0.9, 0.9),
    line.col='black', tube.lty='solid', line.lty='solid', add=TRUE)
{
    if (!missing(df)) {
        pos <- df$pos
        gp.mu <- df$gp.mu
        gp.sd <- df$gp.sd
    }

    # maybe the analytical solution exists, but why derive it
    if (!logspace) {
        cat("transforming to AF space...\n")
        sd.upper <- 1 / (1 + exp(-(gp.mu + sd.mult*gp.sd)))
        sd.lower <- 1 / (1 + exp(-(gp.mu - sd.mult*gp.sd)))
        gp.mu <- 1 / (1 + exp(-gp.mu))
    } else {
        sd.upper <- gp.mu + sd.mult*gp.sd
        sd.lower <- gp.mu - sd.mult*gp.sd
    }

    if (!add) {
        plot(NA, NA, xlim=range(pos), ylim=0:1)
    }

    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)),
        col=tube.col, border=line.col, lty=tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)),
        col=tube.col, border=line.col, lty=tube.lty)
    lines(pos, gp.mu, lwd=2, col=line.col, lty=line.lty)
}
