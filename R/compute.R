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
        dps <- 0:dp
        mut=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd)
        pre1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=2)
        pre2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=2)
        pre <- (pre1 + pre2)/2
        mda1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=4)
        mda2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=4)
        mda <- (mda1 + mda2)/2

        # Step 2: compute model p-values (=alpha) and power to
        # differentiate artifacts from true mutations (=beta) for the
        # two artifact models.
        abc.pv <- sum(mut[mut <= mut[altreads + 1]])
        lysis.pv <- sum(pre[pre <= pre[altreads + 1]])
        lysis.beta <- sum(mut[pre <= pre[altreads + 1]])
        mda.pv <- sum(mda[mda <= mda[altreads + 1]])
        mda.beta <- sum(mut[mda <= mda[altreads + 1]])

        c(abc.pv, lysis.pv, lysis.beta, mda.pv, mda.beta)
    }, altreads, dp, gp.mu, gp.sd)

    rownames(pab) <- c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')

    as.data.frame(t(pab))
}


