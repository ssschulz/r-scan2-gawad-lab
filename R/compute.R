# probability distribution of seeing y variant reads at a site with
# depth d with the estimate (gp.mu, gp.sd) of the local AB. model is
#     Y|p ~ Bin(d, 1/(1+exp(-b))),    b ~ N(gp.mu, gp.sd)
# the marginal distribution of Y (# variant reads) requires integrating
# over b, which has no closed form solution. we approximate it using
# gauss-hermite quadrature. increasing the number of nodes in ghd
# increases accuracy.
# 'factor' allows changing the relationship of ab -> af
# NOTE: for many calls, is best for the caller to compute ghd once
# and supply it to successive calls.

# IMPORTANT!! recomputing this in the function increases runtime
# by approximately 5-fold! Setting this as a global constant is critical!
# XXX: ..but it would be nice if users could change it.
ghd = fastGHQuad::gaussHermiteData(128)
dreads <- function(ys, d, gp.mu, gp.sd, factor=1) { #, ghd=gaussHermiteData(128)) {
    sapply(ys, function(y)
        ghQuad(function(x) {
                # ghQuad requires the weight function on x to be exp(-x^2)
                # (which is NOT a standard normal)
                b <- sqrt(2)*gp.sd*x + gp.mu
                exp(dbinom(y, size=d, prob=1/(factor*(1+exp(-b))), log=TRUE) - log(pi)/2)
            }, ghd
        )
    )
}

mut.model.tables <- function(dp, gp.mu, gp.sd) {
    dps=0:dp
    mut=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd)
    pre1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=2)
    pre2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=2)
    pre <- (pre1 + pre2)/2
    mda1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=4)
    mda2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=4)
    mda <- (mda1 + mda2)/2
    data.frame(dps=dps, mut=mut, pre=pre, mda=mda)
}

# The alpha cutoff:
# For each candidate SNV, we can compute the probability of a
# more extreme event (in terms of the number of alt reads)
# under each artifact model. This roughly corresponds to alpha,
# the type I error rate (false positive rate).
#
# How to determine an alpha cutoff for calling:
# To determine a good calling threshold, we want to know how
# various p-value cutoffs map to false discovery rates (FDRs).
# For example, if our set of candidate SNVs is 99% artifacts,
# then an alpha cutoff of 0.05 would produce an FDR of at
# least 5:1 (~5% of all candidates would be called, 99% of
# which are artifacts; and we don't know how many of the
# remaining true mutations would be called, but it must be
# less than 1% of all candidates). The same 0.05 cutoff would
# perform very differently on a set of candidate mutations
# that is only 1% artifact.
#
# Estimating the FDR:
# The FDR estimate is:
#    FDR = alpha N_A / (alpha N_A + beta N_T).
# N_T and N_A are estimates for the relative rates of true
# mutations and artifacts, respectively, and were calculated
# by comparing VAF distributions between hSNPs and SNV
# candidates. 'alpha' and 'beta' are computed here:
# alpha is the probability of incorrectly calling
# an artifact as a mutation and beta is the probability of
# correctly calling a true mutation.
#
# N.B.:
# Since these distributions are not unimodal, we define "more
# significant" as any event with lower probability.
compute.pvs.and.betas <- function(altreads, dp, gp.mu, gp.sd) {
    pab <- mapply(function(altreads, dp, gp.mu, gp.sd) {
        # Step1: compute dreads for all relevant models:
        # These dreads() calls are the most expensive part of genotyping
        tb <- mut.model.tables(dp, gp.mu, gp.sd)

        # Step 2: compute model p-values (=alpha) and power to
        # differentiate artifacts from true mutations (=beta) for the
        # two artifact models.
        abc.pv <- sum(tb$mut[tb$mut <= tb$mut[altreads + 1]])
        lysis.pv <- sum(tb$pre[tb$pre <= tb$pre[altreads + 1]])
        lysis.beta <- sum(tb$mut[tb$pre <= tb$pre[altreads + 1]])
        mda.pv <- sum(tb$mda[tb$mda <= tb$mda[altreads + 1]])
        mda.beta <- sum(tb$mut[tb$mda <= tb$mda[altreads + 1]])

        c(abc.pv, lysis.pv, lysis.beta, mda.pv, mda.beta)
    }, altreads, dp, gp.mu, gp.sd)

    rownames(pab) <- c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')

    as.data.frame(t(pab))
}


bin.afs <- function(afs, bins=20) {
    sq <- seq(0, 1, 1/bins)
    x <- findInterval(afs, sq, left.open=TRUE)
    # af=0 will be assigned to bin 0 because intervals are (.,.]
    x[x==0] <- 1
    tx <- tabulate(x, nbins=bins)
    names(tx) <- apply(cbind(sq[-length(sq)], sq[-1]), 1, mean)
    tx
}

# the "rough" interval is not a confidence interval. it is just a
# heuristic that demonstrates the range of  reasonably consistent
# somatic mutation burdens
# WARNING! this is the *callable somatic burden*. if you want an
# estimate of the total somatic burden genome-wide, there is
# another function for that provides a better estimate!
estimate.somatic.burden <- function(fc, min.s=1, max.s=5000, n.subpops=10, display=FALSE, rough.interval=0.99) {
    sim <- function(n.muts, g, s, n.samples=1000, diagnose=FALSE) {
        samples <- rmultinom(n=n.samples, size=n.muts, prob=g)
        if (diagnose) {
            boxplot(t(samples))
            lines(s, lwd=2, col=2)
        }
        mean(apply(samples, 2, function(col) all(col <= s)))
    }

    # determine how often a sample of N somatic mutations from the
    # proper het distribution (germlines) "fits" inside the somatic
    # candidate distribution.
    # always try nmut=1-100
    srange <- c(1:100, seq(101, max.s, length.out=n.subpops))
    srange <- srange[srange < max.s]
    fraction.embedded <- sapply(srange, sim, g=fc$g, s=fc$s)
    # max(c(0,... in some cases, there will be no absolutely no overlap
    # between hSNP VAFs and somatic VAFs, leading to fraction.embedded=0
    # for all N.  In this case, return 0 and let the caller decide what
    # to do.
    min.burden <- max(c(0, srange[fraction.embedded >= 1 - (1 -
        rough.interval)/2]))
    max.burden <- max(c(0, srange[fraction.embedded >= (1 - rough.interval)/2]))
    c(min=min.burden, max=max.burden)
}

# {germ,som}.df need only have columns named dp and af
# estimate the population component of FDR
# ignore.100 - ignore variants with VAF=100
fcontrol <- function(germ.df, som.df, bins=20, rough.interval=0.99) {
    germ.afs <- germ.df$af[!is.na(germ.df$af) & germ.df$af > 0]
    som.afs <- som.df$af[!is.na(som.df$af)]
    g <- bin.afs(germ.afs, bins=bins)  # counts, not probabilities
    s <- bin.afs(som.afs, bins=bins)

    # when fcontrol is used on small candidate sets (i.e., when controlling
    # for depth), there may be several 0 bins in s.
    #    g <- g*1*(s > 0)
    # XXX: i don't like this. it doesn't quite match the principle, but
    # rather addresses a limitation of the heuristic method.

    if (length(s) == 0 | all(g == 0))
        return(list(est.somatic.burden=c(0, 0),
             binmids=as.numeric(names(g)),
             g=g, s=s, pops=NULL))

    # returns (lower, upper) bound estimates
    approx.ns <- estimate.somatic.burden(fc=list(g=g, s=s),
        min.s=1, max.s=sum(s), n.subpops=min(sum(s), 100),
        rough.interval=rough.interval)

    cat(sprintf("fcontrol: dp=%d, max.s=%d (%d), n.subpops=%d, min=%d, max=%d\n",
    germ.df$dp[1], nrow(som.df), sum(s), min(nrow(som.df),100),
    as.integer(approx.ns[1]), as.integer(approx.ns[2])))
    pops <- lapply(approx.ns, function(n) {
        #        nt <- pmax(n*(g/sum(g))*1*(s > 0), 0.1)
        nt <- pmax(n*(g/sum(g)), 0.1)
        # ensure na > 0, since FDR would be 0 for any alpha for na=0
        # XXX: the value 0.1 is totally arbitrary and might need to be
        # more carefully thought out.
        na <- pmax(s - nt, 0.1)
        cbind(nt=nt, na=na)
    })

    return(list(est.somatic.burden=approx.ns,
         binmids=as.numeric(names(g)),
         g=g, s=s, pops=pops))
}


estimate.fdr.priors <- function(candidates, hsnps, bins=20, random.seed=0)
{
    # fcontrol -> estimate.somatic.burden relies on simulations to
    # estimate the artifact:mutation ratio.
    cat(sprintf("estimating bounds on number of true mutations in candidate set (seed=%d)..\n",
        random.seed))
    set.seed(random.seed)

    # split candidates by depth; collapse all depths beyond the 80th
    # percentile into one bin
    #hsnps$af <- hsnps$hap1/ hsnps$dp
    max.dp <- as.integer(quantile(hsnps$dp, prob=0.8))
    fcs <- lapply(0:max.dp, function(dp)
        fcontrol(germ.df=hsnps[hsnps$dp == dp,],
                som.df=candidates[candidates$dp == dp,],
                bins=bins)
    )
    fc.max <- fcontrol(germ.df=hsnps[hsnps$dp > max.dp,],
                som.df=candidates[candidates$dp > max.dp,],
                bins=bins)
    fcs <- c(fcs, list(fc.max))

    cat(sprintf("        profiled hSNP and candidate VAFs at depths %d .. %d\n",
        0, max.dp))

    burden <- as.integer(
        c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])),  # min est
          sum(sapply(fcs, function(fc) fc$est.somatic.burden[2])))  # max est
    )

    cat(sprintf("        estimated callable mutation burden range (%d, %d)\n",
        burden[1], burden[2]))
    cat("          -> using MAXIMUM burden\n")

    cat("        estimating true (N_T) and artifact (N_A) counts in candidate set..\n")
    popbin <- ceiling(candidates$af * bins)
    popbin[candidates$dp == 0 | popbin == 0] <- 1

    nt.na <- pbapply::pbmapply(function(dp, popbin) {
        idx = min(dp, max.dp+1) + 1
        if (is.null(fcs[[idx]]$pops))
            c(0.1, 0.1)
        else
            fcs[[idx]]$pops$max[popbin,]
    }, candidates$dp, popbin)

    data.frame(nt=nt.na[1,], na=nt.na[2,])
}




# Finds the (alpha, beta, FDR) that minimizes FDR
min.fdr <- function(pv, alphas, betas, nt, na) {
    x <- data.frame(alpha=alphas, beta=betas, 
        fdr=ifelse(alphas*na + betas*nt > 0,
            alphas*na / (alphas*na + betas*nt), 0))
    x <- rbind(c(1, 0, 1), x)
    x <- x[pv <= x$alpha,] # passing values
    x <- x[x$fdr == min(x$fdr),,drop=F]
    x$fdr[1]
}


# Slower than the current method because the model tables need to
# be fully calculated to find the min. FDR.
# Unlike the real legacy code, the alphas and betas corresponding to
# the minimum FDR are not reported.
compute.fdr.legacy <- function(altreads, dp, gp.mu, gp.sd, nt, na, verbose=TRUE) {
    fdrs <- pbapply::pbmapply(function(altreads, dp, gp.mu, gp.sd, nt, na) {
        # Step1: compute dreads for all relevant models:
        # These dreads() calls are the most expensive part of genotyping
        tb <- mut.model.tables(dp, gp.mu, gp.sd)
    
        # compute BEFORE sorting, so that altreads+1 is the corret row
        lysis.pv <- sum(tb$pre[tb$pre <= tb$pre[altreads + 1]])
        mda.pv <- sum(tb$mda[tb$mda <= tb$mda[altreads + 1]])

        tb <- tb[order(tb$pre),]
        lysis.fdr <- min.fdr(pv=lysis.pv,
            alphas=cumsum(tb$pre), betas=cumsum(tb$mut), nt=nt, na=na)
        
        tb <- tb[order(tb$mda),]
        mda.fdr <- min.fdr(pv=mda.pv,
            alphas=cumsum(tb$mda), betas=cumsum(tb$mut), nt=nt, na=na)

        c(lysis.fdr=lysis.fdr, mda.fdr=mda.fdr)
    }, altreads, dp, gp.mu, gp.sd, nt, na)

    data.frame(lysis.fdr=fdrs[1,], mda.fdr=fdrs[2,])
}


compute.fdr.new <- function(mut.models, fdr.priors) {
    # avoid division by 0
    denom <- mut.models$lysis.pv*fdr.priors$na + mut.models$lysis.beta*fdr.priors$nt
    lysis.fdr <- ifelse(denom > 0, mut.models$lysis.pv*fdr.priors$na / denom, 0)

    denom <- mut.models$mda.pv*fdr.priors$na + mut.models$mda.beta*fdr.priors$nt
    mda.fdr <- ifelse(denom > 0, mut.models$mda.pv*fdr.priors$na / denom, 0)

    data.frame(lysis.fdr=lysis.fdr, mda.fdr=mda.fdr)
}


resample.hsnps <- function(sites, hsnps, M=50, seed=0) {
    # XXX: Random position sampling: maybe add an option to use random
    # instead of somatic candidates?
    # Random positioning is not as realistic as all non-ref sites because
    # it doesn't account for alignability/mappability in the same way as
    # non-ref sites.
    #random.pos <- find.nearest.germline(som=data.frame(chr='X',
    #       pos=sort(as.integer(runif(n=1e5, min=min(spos$pos), max=max(spos$pos))))), 
    #   germ=data, chrs='X')
    #random.pos <- random.pos[abs(random.pos$pos-random.pos$nearest.het)>0,]

    # Distribution of candidates
    # Remember that the nearest hSNP must be on the same chromosome
    tmpsom <- find.nearest.germline(som=sites[order(sites$pos),],
        germ=hsnps,
        chrs=unique(sites$chr))
    spos <- log10(abs(tmpsom$pos-tmpsom$nearest.het))

    # Distribution of hSNP distances
    # Since the training data are already sorted, diff() gives the distance
    # between this SNP and the next. Nearest SNP is the min of the distance
    # to the left and right.
    # Must split by chrom to prevent negative distance between the final
    # site on chr N and the first site on chr N+1
    hsnps$nearest.hsnp <- do.call(c, lapply(unique(hsnps$chr), function(chrom) {
        hsnps <- hsnps[hsnps$chr == chrom,]
        pmin(diff(c(0,hsnps$pos)), diff(c(hsnps$pos, max(hsnps$pos)+1)))
    }))
    hpos <- log10(hsnps$nearest.hsnp)

    # Approximate the distributions
    dist.s <- get.distance.distn(spos)
    dist.h <- get.distance.distn(hpos)
    ds <- dist.s$density[findInterval(hpos, dist.s$breaks, all.inside=T)]
    dh <- dist.h$density[findInterval(hpos, dist.h$breaks, all.inside=T)]
    
    set.seed(seed)
    u <- runif(n=length(ds))
    list(selection=data.frame(dist=hsnps$nearest.hsnp, ds=ds, dh=dh, u=u, keep=u < ds / (M*dh)),
        dist.s=dist.s, dist.h=dist.h)
}
