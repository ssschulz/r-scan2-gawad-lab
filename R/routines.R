# change the pass flag to match a new target FDR
rescore <- function(df, target.fdr, use.pon=FALSE, min.pon.dp=10, quiet=TRUE) {
    newpass <- df$hard.filter &
        df$lysis.fdr <= target.fdr & df$mda.fdr <= target.fdr
    if (use.pon)
        newpass <- newpass & (df$dp >= min.pon.dp & (df$unique.donors <= 1 | df$max.out <= 2))
    if (!quiet)
        cat(sprintf("rescore: %d passing -> %d passing\n", sum(df$pass), sum(newpass)))
    df$pass <- newpass
    df
}


# Read in number of callable basepairs per sample
get.callable <- function(ss.dir, verbose=TRUE) {
    require(yaml)
    ss.config <- file.path(ss.dir, "scan.yaml")
    if (!file.exists(ss.config))
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n",
            ss.config))

    yaml <- read_yaml(ss.config)
    sc.samples <- names(yaml$sc_bams)

    sapply(sc.samples, function(sn) {
        f <- sprintf('%s/callable_regions/%s/callable_regions.bed', ss.dir, sn)
        if (verbose)
            print(f)
        bed <- read.table(f, sep='\t', header=F)
        sum(as.numeric(bed[,3]-bed[,2]))
    })
}



# Read in all genotype data frames
get.scansnv <- function(ss.dir, type='somatic', muttype='snv', verbose=TRUE) {
    if (!(type %in% c('somatic', 'mosaic', 'hsnp_spikein')))
        stop(sprintf("type must be either somatic, mosaic or hsnp_spikein, not '%s'", type))

    if (!(muttype %in% c('snv', 'indel')))
        stop(spritnf("muttype must be either 'snv' or 'indel', not %s", muttype))

    require(yaml)
    ss.config <- file.path(ss.dir, "scan.yaml")
    if (!file.exists(ss.config))
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n",
            ss.config))

    yaml <- read_yaml(ss.config)
    sc.samples <- names(yaml$sc_bams)
    min.sc.alt <- yaml$min_sc_alt

    ret <- lapply(sc.samples, function(s) {
        path.fmt <- "%s_genotypes.rda"
        if (muttype == 'indel' & type == 'somatic')
            path.fmt <- "%s_genotypes.pon_filter.rda"
        f <- file.path(ss.dir, muttype, s, sprintf(path.fmt, type))
        if (verbose)
            print(f)
        load(f)
        # we assume the loaded variable is called 'somatic' below
        # but the mosaic results are called 'mosaic'.
        # just stick it in a variable named somatic anyway. it doesn't matter.
        somatic <- get(ifelse(type == 'hsnp_spikein', 'spikeins', type))

        scalt <- which(colnames(somatic) == make.names(s))+2
        somatic$id <- paste(somatic$chr, somatic$pos, somatic$refnt, somatic$altnt)
        somatic
    })
    names(ret) <- sc.samples

    ret 
}



# given two data frames of somatic and germline locations, annotate
# the somatic data frame with the position of the nearest germline entry.
find.nearest.germline <- function (som, germ, chrs = c(1:22, "X")) {
    som$nearest.het <- NA
    for (chr in chrs) {
        gpos <- germ$pos[germ$chr == chr]
        spos <- som$pos[som$chr == chr]
        gidx <- findInterval(spos, gpos)
        gidx[gidx==0] <- 1  # somatic is to the left of lowest germline, so
                            # the lowest germline is the nearest
        nearest.idx <- ifelse(abs(gpos[gidx] - spos) <= abs(gpos[gidx + 
            1] - spos), gidx, gidx + 1)
        som$nearest.het[som$chr == chr] <- gpos[nearest.idx]
    }
    som
}

# d is a vector of distances to the nearest hSNP
get.distance.distn <- function(d, min=1, max=5) {
    h <- hist(d[d >= min & d <= max], breaks=50, plot=FALSE)
    h$density <- h$density / sum(h$density)
    h
}

# som is the 'somatic' dataframe input to scansnv
# hsnp is the 'data' dataframe from training.rda containing training hSNP sites
# XXX: M=50 is unlikely to be a good value in general
resample.hsnps <- function(som, hsnps, chrom, M=50) {
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
    # 1. only consider candidates for this sample
    #tmpsom <- som[!is.na(som$af) & som$af>0,]
    tmpsom <- find.nearest.germline(som=som[order(som$pos),], germ=hsnps,
        chrs=chrom)
    spos <- log10(abs(tmpsom$pos-tmpsom$nearest.het))

    # Distribution of hSNP distances
    # Since the training data are already sorted, diff() gives the distance
    # between this SNP and the next. Nearest SNP is the min of the distance
    # to the left and right.
    hsnps$nearest.hsnp <- pmin(diff(c(0,hsnps$pos)),
                               diff(c(hsnps$pos, max(hsnps$pos)+1)))
    hpos <- log10(hsnps$nearest.hsnp)

    # Approximate the distributions
    dist.s <- get.distance.distn(spos)
    dist.h <- get.distance.distn(hpos)
    ds <- dist.s$density[findInterval(hpos, dist.s$breaks, all.inside=T)]
    dh <- dist.h$density[findInterval(hpos, dist.h$breaks, all.inside=T)]
    
    u <- runif(n=length(ds))
    list(selection=data.frame(dist=hsnps$nearest.hsnp, ds=ds, dh=dh, u=u, keep=u < ds / (M*dh)),
        dist.s=dist.s, dist.h=dist.h)
}



muttype.map <- c(
    'A>C'='T>G',
    'A>G'='T>C',
    'A>T'='T>A',
    'C>A'='C>A',
    'C>G'='C>G',
    'C>T'='C>T',
    'G>A'='C>T',
    'G>C'='C>G',
    'G>T'='C>A',
    'T>A'='T>A',
    'T>C'='T>C',
    'T>G'='T>G'
)

get.3mer <- function(df) {
    require(BSgenome)
    #require(BSgenome.Hsapiens.UCSC.hg19)
    # we use hs37d5; mostly doesn't matter for autosomes, but chrMT is
    # significantly updated.
    require(BSgenome.Hsapiens.1000genomes.hs37d5)

    comp <- c('A', 'C', 'G', 'T')
    names(comp) <- c('T', 'G', 'C', 'A')

    x <- df

    if (!('muttype' %in% colnames(x))) {
        cat("adding mutation types..\n")
        x$muttype <- muttype.map[paste(x$refnt, x$altnt, sep = ">")]
    }

    #x$ctx <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                    #names=paste("chr", x$chr, sep=''),
    x$ctx <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5,
                    names=x$chr,
                    start=x$pos-1, end=x$pos+1, as.character=TRUE)
    x$ctx.rc <- sapply(strsplit(x$ctx, ""),
                    function(s) paste0(comp[s[c(3,2,1)]], collapse=''))

    x$type.and.ctx <- ifelse(x$refnt == 'C' | x$refnt == 'T',
                       paste0(x$ctx, ":", x$muttype),
                       paste0(x$ctx.rc, ":", x$muttype))
    x
}

# somatic and hsnps must have 'af' and 'dp' columns
get.fdr.tuning.parameters <- function(somatic, hsnps, bins=20, random.seed=0)
{
    cat(sprintf("estimating bounds on somatic mutation rate (seed=%d)..\n",
        random.seed))
    # fcontrol -> estimate.somatic.burden relies on simulations to
    # estimate the artifact:mutation ratio.
    set.seed(random.seed)

    max.dp <- as.integer(quantile(hsnps$dp, prob=0.8))
    fcs <- lapply(0:max.dp, function(dp)
        fcontrol(germ.df=hsnps[hsnps$dp == dp,],
                som.df=somatic[somatic$dp == dp,],
                bins=bins)
    )
    fc.max <- fcontrol(germ.df=hsnps[hsnps$dp > max.dp,],
                som.df=somatic[somatic$dp > max.dp,],
                bins=bins)
    fcs <- c(fcs, list(fc.max))
    cat(sprintf("        profiled hSNP and somatic VAFs at depths %d .. %d\n",
        0, max.dp))

    burden <- as.integer(
        c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])),  # min est
          sum(sapply(fcs, function(fc) fc$est.somatic.burden[2])))  # max est
    )

    cat(sprintf("        estimated callable somatic mutation burden range (%d, %d)\n",
        burden[1], burden[2]))
    cat("          -> using MAXIMUM burden\n")
    list(bins=bins, burden=burden, fcs=fcs, max.dp=max.dp)
}


apply.fdr.tuning.parameters <- function(somatic, fdr.tuning) {
    somatic$popbin <- ceiling(somatic$af * fdr.tuning$bins)
    somatic$popbin[somatic$dp == 0 | somatic$popbin == 0] <- 1

    nt.na <- mapply(function(dp, popbin) {
        idx = min(dp, fdr.tuning$max.dp+1) + 1
        if (is.null(fdr.tuning$fcs[[idx]]$pops))
            c(0.1, 0.1)
        else
            fdr.tuning$fcs[[idx]]$pops$max[popbin,]
    }, somatic$dp, somatic$popbin)

    rownames(nt.na) <- c('nt', 'na')
    nt.na
}

genotype.somatic <- function(gatk, gatk.lowmq, sc.idx, bulk.idx,
    sites.with.ab, sc.cigars, bulk.cigars, cigar.training, cigar.emp.score,
    fdr.tuning, spikein=FALSE,
    cap.alpha=TRUE, cg.id.q=0.05, cg.hs.q=0.05, random.seed=0, target.fdr=0.1,
    bulkref=bulk.idx+1, bulkalt=bulk.idx+2, scref=sc.idx+1, scalt=sc.idx+2,
    min.sc.alt=0, min.sc.dp=0, min.bulk.dp=0)
{
    call.fingerprint <- as.list(environment())
    # almost all of this information is saved in the results
    dont.save <- c('gatk', 'gatk.lowmq', 'sites.with.ab',
        'sc.cigars', 'bulk.cigars')
    call.fingerprint <- call.fingerprint[!(names(call.fingerprint) %in% dont.save)]

    # N.B.: there used to be some randomness in this function, but I
    # believe all of it now resides elsewhere. To be safe, I've left
    # this set.seed so that results can be reproduced.
    set.seed(random.seed)

    cat("step 1: preparing data\n")
    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
    gatk$dp <- gatk[,scalt] + gatk[,scref]
    gatk$af <- gatk[,scalt] / gatk$dp
    gatk$bulk.dp <- gatk[,bulkalt] + gatk[,bulkref]

    # sites only has columns 'chr','pos','refnt','altnt', which match gatk.
    # so this merge call is really just subsetting gatk.
    somatic <- merge(gatk, sites.with.ab, all.y=TRUE)
    cat(sprintf("        %d somatic SNV candidates\n", nrow(somatic)))
    # choose the AB nearest to the AF of each candidate
    # af can be NA if the site has 0 depth
    somatic$gp.mu <- ifelse(!is.na(somatic$af) & somatic$af < 1/2,
        -abs(somatic$gp.mu), abs(somatic$gp.mu))
    somatic$ab <- 1/(1+exp(-somatic$gp.mu))


    cat("step 2: computing p-values for filters\n")
    cat("        allele balance consistency\n")
    somatic$abc.pv <-
        mapply(abc2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, factor=1, dp=somatic$dp)

    cat("        lysis artifacts\n")
    somatic$lysis.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, dp=somatic$dp, div=2)

    cat("        MDA artifacts\n")
    somatic$mda.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu,
            gp.sd=somatic$gp.sd, dp=somatic$dp, div=4)


    cat(sprintf("step 3: tuning FDR = %0.3f\n", target.fdr))
    if (cap.alpha)
        cat(sprintf("        cap.alpha=TRUE: alpha <= %0.3f enforced despite artifact prevalence\n", target.fdr))

    nt.na <- apply.fdr.tuning.parameters(somatic, fdr.tuning)
    somatic <- cbind(somatic, t(nt.na))

    cat("        lysis artifact FDR\n")
    somatic <- cbind(somatic,
                lysis=t(mapply(match.fdr3,
                            pv=somatic$lysis.pv,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=somatic$nt, na=somatic$na,
                            div=2)))

    cat("        MDA artifact FDR\n")
    somatic <- cbind(somatic,
                mda=t(mapply(match.fdr3,
                            pv=somatic$mda.pv,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=somatic$nt, na=somatic$na,
                            div=4)))


    cat("step 5: applying optional alignment filters\n")
    if (missing(gatk.lowmq)) {
        cat("        WARNING: skipping low MQ filters. will increase FP rate\n")
        somatic$lowmq.test <- TRUE
    } else {
        cat(sprintf("        attaching low MQ data..\n"))
        lmq <- gatk.lowmq[,c('chr', 'pos', 'refnt', 'altnt',
                             colnames(gatk.lowmq)[c(scref,scalt,bulkref,bulkalt)])]
        somatic <- merge(somatic, lmq, by=c('chr', 'pos', 'refnt', 'altnt'),
            all.x=TRUE, suffixes=c('', '.lowmq'))
        # At the moment only removing sites that have bulk support when
        # the MQ threshold is lowered.
        cn <- paste0(colnames(gatk.lowmq)[bulkalt], '.lowmq')
        # In spikein mode, known hSNPs are being treated like somatics,
        # so it is necessary to short circuit this test.
        somatic$lowmq.test <- is.na(somatic[,cn]) | somatic[,cn] == 0 | spikein
    }

    somatic <- merge(somatic, sc.cigars, by=c('chr', 'pos'), all.x=T)
    somatic <- merge(somatic, bulk.cigars, by=c('chr', 'pos'), all.x=T,
        suffixes=c('', '.bulk'))
    somatic$id.score.y <- somatic$ID.cigars / somatic$dp.cigars
    somatic$id.score.x <- somatic$ID.cigars.bulk / somatic$dp.cigars.bulk
    somatic$id.score <- cigar.emp.score(training=cigar.training, test=somatic, which='id')
    somatic$hs.score.y <- somatic$HS.cigars / somatic$dp.cigars
    somatic$hs.score.x <- somatic$HS.cigars.bulk / somatic$dp.cigars.bulk
    somatic$hs.score <- cigar.emp.score(training=cigar.training, test=somatic, which='hs')

    cat("        Excessive indel CIGAR ops\n")
    somatic$cigar.id.test <-
        somatic$id.score > quantile(cigar.training$id.score, prob=cg.id.q, na.rm=T)
    cat("        Excessive clipped read CIGAR ops\n")
    somatic$cigar.hs.test <-
        somatic$hs.score > quantile(cigar.training$hs.score, prob=cg.hs.q, na.rm=T)
    

    somatic$dp.test <- somatic$dp >= min.sc.dp & somatic$bulk.dp >= min.bulk.dp

    cat("step 6: calling somatic SNVs\n")
    somatic$hard.filter <- somatic$abc.pv > 0.05 &
        somatic$cigar.id.test & somatic$cigar.hs.test &
        somatic$lowmq.test & somatic$dp.test & somatic[,scalt] >= min.sc.alt
    somatic$pass <- somatic$hard.filter &
        somatic$lysis.fdr <= target.fdr & somatic$mda.fdr <= target.fdr
    cat(sprintf("        %d passing somatic SNVs\n", sum(somatic$pass)))
    cat(sprintf("        %d filtered somatic SNVs\n", sum(!somatic$pass)))

    # return non-data parameters and the output
    # the final RDA files are already large enough to be painful to work with
    return(c(call.fingerprint, list(somatic=somatic)))
}

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
dreads <- function(ys, d, gp.mu, gp.sd, factor=1, ghd=gaussHermiteData(128)) {
    require(fastGHQuad)
    
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

# determine the beta (power) as a function of alpha (FP rate)
# note that this is subject to heavy integer effects, specifically
# when D is low or when af is close to 0 or 1
# div controls the error model
#
# Now generates all possible (alpha,beta) pairs by iterating over all
# possible read supports. No longer needs a "requested alpha" input.
estimate.alphabeta3 <- function(gp.mu, gp.sd, dp=30, div=2) {
    td <- data.frame(dp=0:dp,
        mut=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd),
        err1=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div),
        err2=dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))

    # (td[,4]+td[,5])/2 corresponds to an even mixture of errors on the p
    # and 1-p alleles.
    # XXX: i really don't know how I feel about this mixture. but use it
    # for now, maybe later come up with a better prior for the two error
    # modes.
    td$err <- (td$err1 + td$err2)/2
    td <- td[order(td$err),]
    td$cumerr <- cumsum(td$err)

    # cumerr is the test's alpha, cumsum(td$mut) is the corresponding power
    return(list(td=td, alphas=td$cumerr, betas=cumsum(td$mut)))
}


# NOT equivalent to test(div=1)
abc2 <- function(altreads, gp.mu, gp.sd, dp, factor=1) {
    pmf <- dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=factor)
    p.value <- sum(pmf[pmf <= pmf[altreads + 1]])
    c(p.value=p.value)  # just names the return
}

test2 <- function(altreads, gp.mu, gp.sd, dp, div) {
    # equal mixture of errors in cis and trans
    err <- (dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div) +
            dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))/2

    p.value <- sum(err[err <= err[altreads + 1]])
    c(p.value=p.value)  # just names the return
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
# estimate of the total somatic burden genome-wide, adjust this
# number by the fraction of genome represented in the somatic
# candidate set.
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

# SUMMARY
# given a target FDR, try to match the FDR as best as we can without
# exceeding it.
#
# DETAILS
# for a given (af, dp), find the relation (alpha, beta) to determine
# what sensitivities and specificities are achievable.  using the
# estimated numbers of true mutations nt and artifacts na in the
# population, calculate the corresponding FDRs (alpha,beta,FDR).
# return (alpha, beta, FDR) such that
# FDR = arg max { FDR : FDR <= target.fdr }
# to break ties (when more than one (alpha,beta) produce the same
# FDR), select the triplet with minimal alpha and maximal beta.
# cap.alpha - because there is often such a strong spike at low VAF
#             amongst somatic candidate sites, the FDR control
#             procedure sometimes results in prior error rates of <1%
#             amongst high VAF candidate sSNVs. In this case, it will
#             often accept variants that are very clearly more consistent
#             with the artifact model(s) than with a true mutation. The
#             cap.alpha parameter helps to address this issue by requiring
#             alpha <= target FDR.
#             To be more clear: when nt >> na, FDR matching can pass
#             variants that have very high alpha (e.g., cases with
#             alpha ~ 0.4, beta ~ 0.7). This puts significant pressure
#             on the accuracy of the nt and na estimates.
#
# Now returns the (alpha, beta, fdr) values for the smallest FDR such
# that pv <= alpha.
# THE ONLY DIFFERENCE between this calculation and the published version
# is that there is no cap.alpha. This only kept (a,b,fdr) triplets such
# that alpha <= target.fdr. When nt >> na, the fdr can be severely
# deflated, allowing sites that are clearly concordant with the error
# mode to pass.
match.fdr3 <-
function (pv, gp.mu, gp.sd, dp, nt, na, div = 2)
{
    alphabeta <- estimate.alphabeta3(gp.mu = gp.mu, gp.sd = gp.sd, 
        dp = dp, div = div)
    x <- data.frame(alpha = alphabeta$alphas, beta = alphabeta$betas, 
        fdr = ifelse(alphabeta$alphas * na + alphabeta$betas * 
            nt > 0, alphabeta$alphas * na/(alphabeta$alphas * 
            na + alphabeta$betas * nt), 0))
    x <- rbind(c(1, 0, 1), x)
    x <- x[pv <= x$alpha,] # passing values
    x <- x[x$fdr == min(x$fdr),,drop=F]
    unlist(x[1,])
}
