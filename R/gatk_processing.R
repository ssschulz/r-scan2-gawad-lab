# Some extra work to make sure we only read in the part of the
# table relevant to these two samples. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
read.gatk.table.2sample <- function(path, sc.sample, bulk.sample, region, quiet=FALSE) {
    if (!quiet) cat("Importing GATK table..\n")

    # Step 1: just get the header and detect the columns corresponding to
    # sc.sample and bulk.sample to avoid reading the full matrix.
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
    col.strings <- strsplit(header, '\t')[[1]]
    tot.cols <- length(col.strings)
    sc.idx <- Inf
    if (!missing(sc.sample))
        sc.idx <- which(col.strings == sc.sample)
    bulk.idx <- Inf
    if (!missing(bulk.sample))
        bulk.idx <- which(col.strings == bulk.sample)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            # First 5 are chr, pos, dbsnp, refnt, altnt
            if (i <= 5) {
                cat(sprintf("    (%d)", i), col.strings[i], '\n')
            } else if (any(i == sc.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
            } else if (any(i == bulk.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
            }
        }
    }
    
    # Step 2: really read the tables in, but only the relevant columns
    cols.to.read <- rep("NULL", tot.cols)
    # First 5 are chr, pos, dbsnp, refnt, altnt
    cols.to.read[1:5] <- c('character', 'integer', rep('character', 3))
    # Read 3 columns for the single cell, 3 columns for bulk
    if (!missing(sc.sample))
        cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
    if (!missing(bulk.sample))
        cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')

    gatk <- read.tabix.data(tf=tf, region=region, header=header, quiet=quiet, colClasses=cols.to.read)
    close(tf)

    if (!missing(sc.sample)) {
        new.sc.idx <- which(colnames(gatk) == sc.sample)
        colnames(gatk)[new.sc.idx+1:2] <- c('scref', 'scalt')
    }

    if (!missing(bulk.sample)) {
        new.bulk.idx <- which(colnames(gatk) == bulk.sample)
        colnames(gatk)[new.bulk.idx+1:2] <- c('bref', 'balt')
    }

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- col.strings[1:5]
    if (!missing(sc.sample))
        cols.to.keep <- c(cols.to.keep, sc.sample, c('scref', 'scalt'))
    if (!missing(bulk.sample))
        cols.to.keep <- c(cols.to.keep, bulk.sample, c('bref', 'balt'))

    gatk <- gatk[,..cols.to.keep]

    gatk
}


# annotate `gatk.meta` with 'bulk.gt', 'bulk.dp', 'bulk.af', 'bref' and 'balt' taken from `gatk`,
# which are the bulk genotype string assigned by HaplotypeCaller, bulk depth, bulk VAF, the number
# of ref supporting bulk reads and mutation supporting bulk reads.
#
# only `gatk.meta` is modified (by reference).
annotate.gatk.bulk <- function(gatk.meta, gatk, bulk.sample, quiet=FALSE) {
    #header <- colnames(gatk)
    #col.strings <- strsplit(header, '\t')[[1]]
    col.strings <- colnames(gatk)
    tot.cols <- length(col.strings)
    bulk.idx <- Inf
    if (!missing(bulk.sample))
        bulk.idx <- which(col.strings == bulk.sample)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (any(i == bulk.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
            }
        }
    }

    gatk.meta[, c('bulk.gt', 'bref', 'balt') :=
        list(gatk[[bulk.idx]], gatk[[bulk.idx+1]], gatk[[bulk.idx+2]])]
    gatk.meta[, c('bulk.dp', 'bulk.af') := bref+balt, balt/(bref+balt)]
}


# `gatk` is a data.table, so all of these updates happen by reference.
# no need to return the result.
#
# This can be slow with add.mutsig, particularly for indels. Best used
# on chunked data.
annotate.gatk <- function(gatk, genome.string, genome.object, add.mutsig=TRUE) {
    data.table::setkey(gatk, chr, pos, refnt, altnt)

    # Determine SNV/indel status and then annotate mutation signature channels
    gatk[, muttype := ifelse(nchar(refnt) == 1 & nchar(altnt) == 1, 'snv', 'indel')]
    # allow for fast selection of SNVs or indels
    setindex(gatk, muttype)

    if (add.mutsig) {
        gatk[muttype == 'snv',
            mutsig := get.3mer(chr=chr, pos=pos, refnt=refnt, altnt=altnt, genome=genome.object)]
        chs <- classify.indels(gatk[muttype == 'indel'], genome.string=genome.string)
        gatk[muttype == 'indel', mutsig := chs]
    }

    # FIXME: legacy SCANSNV behavior did indeed use every alt read count
    # other than the named bulk sample to determine somatic candidate sites.
    # This does not account for jointly analyzing other bulk samples (that
    # perhaps were not the primary bulk to be compared against) or other
    # single cells that were not included in analysis for whatever reason
    # (like low-depth cells).
    # We are continuing legacy behavior, but I'd prefer to fix it at some
    # point.
    alts <- which(colnames(gatk) == 'alt')
    sum.alts <- rowSums(as.matrix(gatk[,..alts])) - gatk$balt

    # N.B. rowSums requred >= min.sc.alt, but in legacy uses this was always 2.
    # FIXME: it'd be good to expose this to the end user. Especially for mosaic
    # mutation calling.
    gatk[, somatic.candidate := balt == 0 & bulk.gt == '0/0' & dbsnp == '.' & sum.alts >= 2]
}


# 'gatk' is a data.table to be modified by reference
annotate.gatk.lowmq <- function(gatk, path, bulk, region, quiet=FALSE) {
    lowmq <- read.gatk.table.2sample(path, bulk.sample=bulk, region=region, quiet=quiet)
    data.table::setkey(lowmq, chr, pos, refnt, altnt)

    gatk[lowmq, on=.(chr,pos,refnt,altnt), balt.lowmq := i.balt]
}


# 'gatk' is a data.table to be modified by reference
# 'phasing.path' points to a phased VCF output by either SHAPEIT or EAGLE
# with sample column 'phasedgt'. The VCF file must have a header line
# beginning with #CHROM, as is usual for the VCF spec.
#
# IMPORTANT: this VCF contains no information about single cells, only bulk.
#
# to allow the very simple table join to work well, we assume phasing was
# restricted to biallelic sites for this individual. while not ideal, the
# performance of such phasing has been good enough.
#
# To annotate phasing status (phased.hap1, phased.hap2) the estimated
# phase needs to be joined to a single cell table with alt and ref read
# counts. Based on phase, (alt,ref) maps to either (hap1,hap2) or (hap2,hap1).
annotate.gatk.phasing <- function(gatk, phasing.path, region, quiet=FALSE) {
    phase.data <- read.tabix.data(path=phasing.path, region=region, quiet=quiet)

    if (ncol(phase.data) != 10)
        stop('expected single-sample VCF format with 10 columns')

    # phasing.path is a single sample standard VCF, so we specify the column
    # format explicitly.
    colnames(phase.data) <- c('chr', 'pos', 'dbsnp', 'refnt', 'altnt', 'qual', 'filter', 'info', 'format', 'phasedgt')

    # This assumes "GT" is the first element of the GT string format, which isn't
    # guaranteed but is the case for our data.

    unrecognized <- setdiff(unique(phase.data$phasedgt), c('1|0', '0|1', './.'))
    if (length(unrecognized) > 0)
        stop(paste('phasing genotypes expected to be either 0|1, 1|0 or ./., but found', unrecognized, collapse='\n'))

    # First join the phase genotype (string is either 0|1, 1|0 or ./. if no call)
    gatk[phase.data, on=.(chr,pos,refnt,altnt), phased.gt := i.phasedgt]
}


rename.and.use.later <- function()
{
    gatk[, c('phased.hap1', 'phased.hap2') :=
        list(ifelse(phased.gt == '0|1', scref, scalt),
             ifelse(phased.gt == '0|1', scalt, scref))]

    if (parsimony.phasing)
        adjust.phase(gatk, dist.cutoff=parsimony.dist.cutoff, quiet=quiet)
}


# Mark the germline hSNPs and hIndels that should be used as training sites.
# Note that hIndels aren't used for AB model fitting; currently they are used
# to build CIGAR op filters, determine the FDR prior distns and for sensitivity
# estimation used in total burden calculations.
annotate.gatk.training <- function(gatk, single.cell, bulk) {
    bulk.gt <- gatk[[bulk]]
    sc.gt <- gatk[[single.cell]]
    gatk[, training.site := phased.gt != './.' & sc.gt != './.' & bulk.gt != './.']
}


# Again, modifying 'gatk' by reference.
# Returns list of resampling auxiliary data with one entry for SNVs and one for indels
# and modifies `gatk` by reference.
# FIXME: germline indels, not germline SNVs, were used by legacy SCAN2 calling when
# downsampling. This is right and wrong in various scenarios:
#    - For sensitivity estimation, this is WRONG because only hSNPs are used to
#      in AB model inference. Thus, to best reflect sensitivity due to proximity
#      to AB model training sites (somatic sites closer to training sites have more
#      accurate AB estimates and therefore probably higher sensitivity), only hSNPs
#      should be used.
#    - For CIGAR op filtering, this is CORRECT because somatic indel CIGAR profiles
#      should be compared to germline indel CIGAR profiles. Of course, applying certain
#      CIGAR op filters (indel ops) is actually incorrect either way.
# In any case, this likely has a relatively insignificant effect so we are leaving it as
# it was in legacy calling for now.
gatk.resample.phased.sites <- function(gatk, M=20, seed=0)
{
    ret <- list()
    for (mt in c('snv', 'indel')) {
        aux.data <- resample.germline(
            sites=gatk[somatic.candidate == TRUE & muttype == mt],
            hsnps=gatk[training.site == TRUE & muttype == mt],
            M=M, seed=seed)

        # aux.data$selection is aligned to the input 'hsnps' table
        gatk[training.site == TRUE & muttype == mt, resampled.training.site := aux.data$selection$keep]
        ret[[mt]] <- aux.data
    }

    gatk[is.na(resampled.training.site), resampled.training.site := FALSE]
    setindex(gatk, resampled.training.site)

    return(ret)
}
