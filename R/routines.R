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
    ss.config <- file.path(ss.dir, "scan.yaml")
    if (!file.exists(ss.config))
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n",
            ss.config))

    yaml <- yaml::read_yaml(ss.config)
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

    ss.config <- file.path(ss.dir, "scan.yaml")
    if (!file.exists(ss.config))
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n",
            ss.config))

    yaml <- yaml::read_yaml(ss.config)
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
resample.hsnps.old <- function(som, hsnps, chrom, M=50) {
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

old.get.3mer <- function (df, genome = BSgenome.Hsapiens.1000genomes.hs37d5) 
{
    if ("type.and.ctx" %in% colnames(df)) 
        return(df)
    require(BSgenome)
    x <- df
    if (!("muttype" %in% colnames(x))) {
        cat("adding mutation types..\n")
        x$muttype <- muttype.map[paste(x$refnt, x$altnt, sep = ">")]
    }
    tmp <- getSeq(genome, GRanges(seqnames = x$chr, ranges = IRanges(start = x$pos - 
        1, end = x$pos + 1)))
    x$ctx <- as.character(tmp)
    x$ctx.rc <- as.character(reverseComplement(tmp))
    x$type.and.ctx <- ifelse(x$refnt == "C" | x$refnt == "T", 
        paste0(x$ctx, ":", x$muttype), paste0(x$ctx.rc, ":", 
            x$muttype))
    x
}


get.3mer <- function(df, chr, pos, refnt, altnt, genome=BSgenome.Hsapiens.1000genomes.hs37d5,verbose=FALSE) {
    if (!missing(df)) {
        chr <- df$chr
        pos <- df$pos
        refnt <- df$refnt
        altnt <- df$altnt
    } else {
        if (missing(chr) | missing(pos) | missing(refnt) | missing(altnt))
            stop('must supply either "df" or all of chr, pos, refnt, alt')
    }
    
    muttype <- muttype.map[paste0(refnt, '>', altnt)]

    tmp <- getSeq(genome, GRanges(seqnames=chr,
                    ranges=IRanges(start=pos-1, end=pos+1)))

    ctx <- as.character(tmp)
    ctx.rc <- as.character(reverseComplement(tmp))
    type.and.ctx <- paste0(ifelse(refnt == 'C' | refnt == 'T', ctx, ctx.rc), ':', muttype)
    # sbs96 converts the string form into a factor to enforce ordering
    sbs96(type.and.ctx)
}


# Extremely simple heuristic for improving short-range phasing
# Ties caused by the abs() > abs() statement produce a bias in
# adjusted AFs.
# 10kb is probably much larger than the correlation signal
# caries, but it doesn't actually matter. 
#
# IMPORTANT: adjust.phase does NOT perform accurate long-range
# phasing. By setting the phase to the most agreeable VAF of
# the nearest hSNP, the phase for small ranges of VAFs near
# 0.5 will likely be incorrect; however, these small phase
# swaps are unlikely to affect the mutation models.
# The intention is to resolve regions with high imbalance and
# to thereby prevent AB "ping-ponging" between ABs near 0 and 1.
#
# N.B. with enough single cells, a method that harmonizes VAFs
# for multiple cells simultaneously would perform accurate long
# range phasing.
#
# 'pht' is a data.table that is modified by reference.
adjust.phase <- function(pht, dist.cutoff=1e4, quiet=FALSE) {
    if (!quiet)
    cat("WARNING: adjust.phase is EXPERIMENTAL!\n")

    pos <- pht$pos
    dp <- pht$phased.hap1 + pht$phased.hap2
    # copy will be overwritten as the sapply below runs
    af <- pht$phased.hap1 / dp
    cat(sprintf("adjust.phase: neighboring hSNP correlation before adjustment: %0.3f\n",
        cor(af[-1], af[-nrow(pht)], use='complete.obs')))

    # this can't be by independently evaluating each site. after evaluation of one
    # site, the next site's results may change.
    swapped <- rep(F, length(af))
    sapply(2:length(af), function(i) {
        swap.i <- pos[i] - pos[i-1] < dist.cutoff &
                !is.na(af[i-1]) & !is.na(af[i]) &
                abs(af[i-1]-af[i]) > abs(af[i-1]-(1-af[i]))
        swapped[i] <<- swap.i
        if (swap.i)
            af[i] <<- 1-af[i]
    })

    pht[, c('phased.hap1', 'phased.hap2', 'phase.adjusted') :=
        list(ifelse(swapped, phased.hap2, phased.hap1),
             ifelse(swapped, phased.hap1, phased.hap2),
             swapped)]

    if (!quiet) {
        cat(sprintf("adjust.phase: neighboring hSNP correlation after adjustment: %0.3f\n",
            cor(af[-1], af[-nrow(pht)], use='complete.obs')))
    }
}
