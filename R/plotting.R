id83.cols <- col <- rep(c('#FBBD75', '#FC7F24', '#B0DB8E', '#3B9F36',
    '#FBC9B6', '#F8896D', '#EE453A', '#B91C22', '#C5D4E4', '#8DBAD2',
    '#4D98C6', '#1D65A8', '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'),
    c(rep(6,12), 1,2,3,5))
id83.channel.order <- paste(
    c(rep(1,24), rep(rep(2:5,each=6),2),c(2,3,3,4,4,4,5,5,5,5,5)),
    c(rep(c('Del', 'Ins'), each=12), rep(c('Del', 'Ins'), each=24), rep('Del', 11)),
    c(rep(rep(c('C','T'), each=6), 2), rep('R',48), rep('M', 11)),
    c(rep(0:5, 12), c(1, 1,2, 1:3, 1:5)),
    sep=':')


sbs96.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey',
    'chartreuse3', 'pink2'), each=16)
sbs96.channel.order <- paste0(rep(c("A", "C", "G", "T"), each = 4), rep(c("C", "T"), 
    each = 48), rep(c("A", "C", "G", "T"), times = 4), ":", rep(c("C", "T"), 
    each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
    rep(c("A", "C", "G"), each = 16)))

sbs96 <- function(x, colname) {
    factorize.mutsig(x, sbs96.channel.order)
}

id83 <- function(x) {
    factorize.mutsig(x, id83.channel.order)
}

factorize.mutsig <- function(x, channel.order) {
    factor(x, levels=channel.order, ordered=TRUE)
}

# Simple extension of R's base table() function to make normalizing
# mutation spectra more convenient.
# 
# x - a factor()ized vector of mutation signature channels
# eps - a minimum count value to add to all channels in the spectrum.
#       added *BEFORE* conversion to fraction, if fraction=TRUE. eps is
#       ignored if fraction=FALSE.
# fraction - express the mutation spectrum as a probability density
#       function rather than counts; i.e., ensure the spectrum sums to 1.
as.spectrum <- function (x, eps = 0.1, fraction = TRUE) {
    t <- table(x)  # when acting on factors, table() will report 0 counts
                   # for channels that have 0 mutations
    if (fraction) {
        t <- t + eps
        t <- t/sum(t)
    }
    t
}


# x - a vector of id83() factors OR a SCAN2 object
# spectrum - an already-tabulated id83 factor spectrum
# either x or spectrum can be supplied, but not both
plot.sbs96 <- function(x, spectrum, xaxt='n', legend=FALSE, ...) {
    if (missing(x) & missing(spectrum))
        stop('exactly one of "x" or "spectrum" must be supplied')

    if (!missing(x) & is(x, 'SCAN2'))
        x <- sbs96(x@gatk$mutsig)   # Indels are automatically ignored because they don't match any of the known SBS96 channels

    if (missing(spectrum))
        spectrum <- table(x)

    p <- barplot(spectrum, las=3, col=sbs96.cols,
        space=0.5, border=NA, xaxt=xaxt, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    if (legend) {
        # mutation types are [context]:[refbase]>[altbase]
        legend('topright', ncol=2, legend=c('C>A','C>G','C>T','T>A','T>C','T>G'),
            fill=sbs96.cols[seq(1, length(sbs96.cols), 16)])
    }
}


# x - a vector of id83() factors OR a SCAN2 object
# spectrum - an already-tabulated id83 factor spectrum
# either x or spectrum can be supplied, but not both
#
# detailed.x.labels - annotate each bar in the barplot with the full
#      mutation class. E.g., "1:Del:C:0". When plotting many separate
#      panels over X11, this can be very slow.
plot.id83 <- function(x, spectrum, proc, xaxt='n',
    col, border, detailed.x.labels=FALSE, ...) {

    if (missing(x) & missing(spectrum))
        stop('exactly one of "x" or "spectrum" must be supplied')

    if (!missing(x) & is(x, 'SCAN2'))
        x <- id83(x@gatk$mutsig)  # SNVs are automatically ignored because they don't match any of the known ID83 channels

    if (missing(spectrum))
        spectrum <- table(x)
    # else it's already tabulated

    if (missing(border))
        border <- id83.cols
    if (missing(col))
        col <- id83.cols

    x.names <- if (detailed.x.labels) names(spectrum) else ''
    p <- barplot(spectrum, las = 3, col = col, names.arg = x.names,
        space = 0.5, border = border, cex.names=0.7, xaxt=xaxt, ...)
    abline(v = (p[seq(6, length(p) - 11, 6)] + p[seq(7, length(p)-10,6)])/2, col="grey")

    if (xaxt != 'n') {
        mtext(text=c('del C', 'del T', 'ins C', 'ins T', 'del 2', 'del 3', 'del 4',
            'del 5+', 'ins 2', 'ins 3', 'ins 4', 'ins 5+', 'microhom.'),
            side=1, at=c(mean(p[1:6]), mean(p[7:12]), mean(p[13:18]),
                mean(p[19:24]), mean(p[25:30]), mean(p[31:36]),
                mean(p[37:42]), mean(p[43:48]), mean(p[49:54]),
                mean(p[55:60]), mean(p[61:66]), mean(p[67:72]), mean(p[73:83])))
    }
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
setGeneric("plot.region", function(object, site=NA, chrom=NA, position=NA, upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=TRUE, show.all.candidates=FALSE)
    standardGeneric("plot.region"))
setMethod("plot.region", "SCAN2", function(object, site=NA, chrom=NA, position=NA,
    upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=TRUE, show.all.candidates=FALSE)
{
    check.slots(object, c('gatk', 'ab.estimates'))
    if (recompute)
        check.slots(object, 'ab.fits')

    if (!missing(site) & (!missing(chrom) | !missing(position)))
        stop("either site or (chrom,pos) must be specified, but not both")
    if ((!missing(chrom) & missing(position)) | (missing(chrom) & !missing(position)))
        stop("both chrom and pos must be specified")

    if (!missing(site)) {
        chrom <- object@gatk[site,]$chr
        position <- object@gatk[site,]$pos
        cat('using site', site, '\n')
        print(object@gatk[site,])
    }

    # Sites at which AB was estimated in the genotyper.
    # Each site in 'd' will be plotted with a point.
    # Don't keep sites with no reads in this sample, unless it's a training
    # site. The majority of these 0 read, non-germline sites are here
    # because they had reads in a different sample.
    d <- object@gatk[chr == chrom &
        pos >= position - upstream & pos <= position + upstream &
        (training.site | somatic.candidate) & bulk.gt != '1/1']

    # Recompute does a better job of showing the model, but it does
    # cost some CPU.
    if (recompute) {
        cat("estimating AB in region with", gp.extend, " basepair flanks..\n")
        # ensure that we estimate at exactly pos and at all sites in 'd'
        # other loci are not sites reported in the GATK table, they are
        # only there to make smooth lines.
        est.at <- c(seq(position - upstream, position-1, length.out=n.gp.points/2), position,
                    seq(position+1, position + downstream, length.out=n.gp.points/2),
                    d$position)
        est.at <- sort(unique(est.at))
        fit.chr <- object@ab.fits[chrom,,drop=FALSE]
        newdt <- object@gatk[training.site == TRUE & chr == chrom]
        # infer.gp requires hsnps to have hap1 and hap2 columns
        newdt[, c('hap1', 'hap2') := list(phased.hap1, phased.hap2)]
        gp <- infer.gp1(ssnvs=data.frame(chr=chrom, pos=est.at),
            fit=fit.chr, hsnps=newdt, flank=gp.extend, max.hsnps=150)
        gp <- data.frame(chr=chrom, pos=est.at, ab=1/(1+exp(-gp[,'gp.mu'])), gp)
    } else {
        # Rely on precomputed GP mu/sd (relevant for ALLSITES mode)
        gp <- d
        # Flip the allele balance to match VAF, as is done when computing
        # the AB true and artifact models.
        gp$gp.mu <- match.ab(af=gp$af, gp.mu=gp$gp.mu)
    }

    plot.gp.confidence(df=gp, add=FALSE,
        xlab=paste('Chrom', d$chr[1], 'position'),
        ylab='Allele fraction')
    points(d[training.site==TRUE]$pos,
        d[training.site==TRUE, phased.hap1/(phased.hap1+phased.hap2)],
        pch=20, cex=1, col=1, ylim=c(-0.2,1))

    # 5*max : restrict the depth to the bottom 20% of plot
    lines(d$pos, d$dp/(5*max(d$dp)), type='h', lwd=2)
    text(d$pos[which.max(d$dp)], 1/5, max(d$dp))
    abline(h=30/(5*max(d$dp)), lty='dotted', col='grey')

    lines(gp$pos, gp$ab/2, lwd=2, col=2)
    lines(gp$pos, (1-gp$ab)/2, lwd=2, col=2)

    # If a site was given, emphasize it
    if (!missing(site)) {
        abline(v=position, lty='dotted')
        abline(h=d[pos==position]$af, lty='dotted')
        points(position, d[pos == position]$af, pch=4, cex=1.5, lwd=2, col=3)
    }

    if (show.all.candidates) {
        legend('topright', legend=c('Training hSNP', 'Candidate', 'Target site'),
            pch=c(20,20,4), col=1:3, pt.cex=c(1,1.5,1.5), bty='n')
        # plot everything except the requested site
        points(d[pos != position & somatic.candidate == TRUE, .(pos, af)], pch=20, col=2)
    } else { 
        legend('topright', legend=c('Training hSNP', 'Target site'),
            pch=c(20,4), col=c(1,3), pt.cex=c(1,1.5))
    }
})

# Add 95% confidence bands for AB model to plot.
# NOTE: the 95% probability interval is in the NORMAL space, not the
# fraction space [0,1]. After the logistic transform, the region in
# fraction space may not be an HDR.
plot.gp.confidence <- function(pos, gp.mu, gp.sd, df, sd.mult=2,
    logspace=FALSE, tube.col=rgb(0.9, 0.9, 0.9),
    line.col='black', tube.lty='solid', line.lty='solid', add=TRUE, ...)
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
        plot(NA, NA, xlim=range(pos), ylim=0:1, ...)
    }

    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)),
        col=tube.col, border=line.col, lty=tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)),
        col=tube.col, border=line.col, lty=tube.lty)
    lines(pos, gp.mu, lwd=2, col=line.col, lty=line.lty)
}


plot.abmodel.covariance <- function(object, bin.breaks=c(1, 10^seq(1,5,length.out=50)), ...) {
    plot.mle.fit <- function(object, ...) {
        ps <- colMeans(object@ab.fits)
        a=ps['a']; b=ps['b']; c=ps['c']; d=ps['d']
        curve(K.func(x, y=0, a=a, b=b, c=c, d=d)/(exp(a)+exp(c)), ...)
    }

    # Use finer binning than the standard call
    neighbor.approx <- approx.abmodel.covariance(object, bin.breaks=bin.breaks)

    # [-1] - use right-hand side of interval for plotting
    plot(neighbor.approx[,.(max.d, corrected.cor)],
        log='x', type='b', pch=16, ylim=0:1,
        xlab='Distance between hSNPs (log10)', ylab='Correlation between hSNP VAFs', ...)
    lines(neighbor.approx[,.(max.d, observed.cor)],
        type='b', pch=1, lty='dotted')
    plot.mle.fit(object, add=TRUE, col=2, lwd=2)
    legend('topright', pch=c(16,1,NA), lwd=c(1,1,2), col=c(1,1,2), lty=c('solid','dotted','solid'),
        legend=c('Adjacent hSNP approx. (corrected)', 'Adjacent hSNP approx. (observed)', 'MLE fit (avg. over chroms)'))
    abline(v=c(150,300), lty='dotted')
}


plot.depth.profile <- function(object, keep.zero=FALSE, quantile=0.99, ...) {
    require(viridisLite)
    # row and column 1 correspond to 0 depth. these usually completely
    # drown out the rest of the depth signal.
    d <- object@depth.profile$dptab
    maxdp <- object@depth.profile$clamp.dp
    x=0:maxdp
    y=0:maxdp
    if (!keep.zero) {
        d <- d[-1,][,-1]
        x <- x[-1]
        y <- y[-1]
    }

    # Cut the plot down to 95% (=quantile option) of the genome in each direction
    xmax <- which(cumsum(rowSums(d))/sum(d) >= quantile)[1]
    if (length(xmax) == 0)  # if not found, take the whole thing
        xmax <- maxdp
    ymax <- which(cumsum(colSums(d))/sum(d) >= quantile)[1]
    if (length(ymax) == 0)  # if not found, take the whole thing
        ymax <- maxdp

    image(x=x, y=y, d, col=viridisLite::viridis(100),
        xlim=c(0,xmax), ylim=c(0,ymax),
        xlab=paste(names(dimnames(d))[1], ' (single cell) depth'),
        ylab=paste(names(dimnames(d))[2], ' (buk) depth'))
}


# tilewidth=1kb by default. the het germline SNP rate in humans is about 1/1.5kb. so
# to get a min. genomic region containing roughly ~100 hSNPs, need 150kb = 150 tiles.
plot.sensitivity <- function(object, min.tiles=150) {
    layout(t(1:2))
    for (mt in c('snv', 'indel')) {
        maj <- assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='maj')
        min <- assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='min')
        plot(maj[n > min.tiles, .(pred, sens)], lwd=2, type='b', pch=20, col=1, main=mt, ylim=0:1,
            xlab="Predicted sensitivity based on local covariates",
            ylab='Actual sensitivity for germline het sites')
        lines(min[n > min.tiles, .(pred, sens)], lwd=2, type='b', pch=20, col=2)
        abline(coef=0:1)
        legend('topleft', lwd=2, col=1:2, legend=c("Major allele", 'Minor allele'))
    }
}


plot.sensitivity.covs <- function(object, muttype=c('snv', 'indel'),
    covs=c('gp.mu', 'gp.sd', 'mean.sc.dp', paste0(muttype, '.', 'n.training.neighborhood'),
           paste0('bases.gt.', muttype, '.sc.min.dp')), min.tiles=150) {
    muttype <- match.arg(muttype)

    layout(matrix(1:(2*length(covs)), nrow=2))
    par(mar=c(5,4,1,1))
    for (cov in covs) {
        xform <- identity
        if (cov == 'gp.mu' | cov == 'gp.sd')  # abs() does nothing for gp.sd, which is >0.
            xform <- function(x) round(abs(x),2)
        if (cov == 'mean.sc.dp')
            xform <- function(x) round(x,0)   # integerize, but nearest

        pred.maj <- paste0('pred.', muttype, '.maj')
        pred.min <- paste0('pred.', muttype, '.min')
        # is.na(mean.bulk.dp) is essentially an alias for tiles in unassembled genome
        # regions. there are ~71 tiles containing hSNPs with is.na(mean.bulk.dp) vs.
        # 195,046 tiles with no hSNPs in the neighborhood of +/- 10kb.
        tab <- object@spatial.sensitivity$data[!is.na(mean.bulk.dp),
            .(sens.maj=mean(get(pred.maj), na.rm=TRUE),
              sens.min=mean(get(pred.min), na.rm=TRUE),
              n.tiles=nrow(.SD)),
            by=.(cov=xform(get(cov)))][order(cov)][n.tiles >= min.tiles]
        plot(tab[,.(cov, sens.maj)], pch=20, ylim=0:1,
            xlab=cov, ylab='Predicted sensitivity')
        lines(tab[,.(cov, sens.maj)])
        points(tab[,.(cov, sens.min)], pch=20, col=2)
        lines(tab[,.(cov, sens.min)], col=2)
        legend('bottomright', legend=c('Major allele', 'Minor allele'), col=1:2, pch=20, lwd=1)
        plot(tab[,.(cov, 100*n.tiles/sum(n.tiles))], pch=20,
            xlab=cov, ylab='Percent of genome')
    }
}
