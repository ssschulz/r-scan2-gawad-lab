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
    sc.idx <- which(col.strings == sc.sample)
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
    cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
    cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')
    gatk <- read.tabix.data(tf=tf, region=region, header=header, quiet=quiet, colClasses=cols.to.read)
    close(tf)
    new.sc.idx <- which(colnames(gatk) == sc.sample)
    new.bulk.idx <- which(colnames(gatk) == bulk.sample)
    colnames(gatk)[new.sc.idx+1:2] <- c('scref', 'scalt')
    colnames(gatk)[new.bulk.idx+1:2] <- c('bref', 'balt')

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- c(col.strings[1:5], sc.sample, c('scref', 'scalt'), bulk.sample, c('bref', 'balt'))
    gatk <- gatk[,..cols.to.keep]

    gatk
}


# `gatk` is a data.table, so all of these updates happen by reference.
# no need to return the result.
#
# This can be slow with add.mutsig, particularly for indels.
annotate.gatk <- function(gatk, add.mutsig=TRUE) {
    # Add some convenient calculations
    # Not sure if it's necessary to split these calculations up.
    gatk[, c('dp', 'bulk.dp') :=
        list(scalt+scref, bref+balt)]
    gatk[, c('af', 'bulk.af') :=
        list(scalt/dp, balt/bulk.dp)]
    data.table::setkey(gatk, chr, pos, refnt, altnt)

    # Determine SNV/indel status and then annotate mutation signature channels
    gatk[, muttype := ifelse(nchar(refnt) == 1 & nchar(altnt) == 1, 'snv', 'indel')]
    # allow for fast selection of SNVs or indels
    setindex(gatk, muttype)

    if (add.mutsig) {
        gatk[muttype == 'snv',
            mutsig := get.3mer(chr=chr, pos=pos, refnt=refnt, altnt=altnt, genome=object@genome.object)]
        chs <- classify.indels(gatk[muttype == 'indel'], genome.string=object@genome.string)
        gatk[muttype == 'indel', mutsig := chs]
    }
}


# 'gatk' is a data.table to be modified by reference
annotate.gatk.lowmq <- function(gatk, path, single.cell, bulk, region, quiet=FALSE) {
    lowmq <- read.gatk.table.2sample(path, single.cell, bulk, region=NULL, quiet=quiet)
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
annotate.gatk.phasing <- function(gatk, phasing.path) {
    # 'skip=str' in fread actually means skip to the first line containing str
    phase.data <- fread(phasing.path, skip='#CHROM', header=TRUE,
        colClasses=c(`#CHROM`='character'))
    colnames(phase.data)[1:5] <- c('chr', 'pos', 'dbsnp', 'refnt', 'altnt')

    # This assumes "GT" is the first element of the GT string format, which isn't
    # guaranteed but is the case for our data.
    phase.data$phasedgt <- sapply(strsplit(phase.data$phasedgt, ':'), head, 1)

    unrecognized <- setdiff(unique(phase.data$phasedgt), c('1|0', '0|1', './.'))
    if (length(unrecognized) > 0)
        stop(paste('phasing genotypes expected to be either 0|1, 1|0 or ./., but found', unrecognized, collapse='\n'))

    # First join the phase genotype (string is either 0|1, 1|0 or ./. if no call)
    gatk[phase.data, on=.(chr,pos,refnt,altnt), phased.gt := i.phgt]

    gatk[, c('phased.hap1', 'phased.hap2') :=
        list(ifelse(phased.gt == '1|0', scalt, scref),
             ifelse(phased.gt == '0|1', scref, scalt))]
}


# Mark the germline hSNPs and hIndels that should be used as training sites.
# Note that hIndels aren't used for AB model fitting; currently they are used
# to build CIGAR op filters, determine the FDR prior distns and for sensitivity
# estimation used in total burden calculations.
annotate.gatk.training <- function(gatk, single.cell, bulk) {
    bulk.gt <- gatk[[bulk]]
    sc.gt <- gatk[[single.cell]]
    gatk[, training.site := phased.gt != './.' & sc.gt != './.' & bulk.gt == './.']
}


# Make sure the supplied GATK table has been fully annotated.
# Intended to be used by the SCAN2 object.
check.gatk <- function(gatk) {
    check.one <- function(cn, action)
        if (!(cn %in% colnames(gatk)))
            stop(paste(cn, 'not detected in gatk, please', action))

    check.one('balt.lowmq', 'join mmq1 data')
    check.one('phased.gt', 'join phasing data')
    check.one('phased.hap1', 'join phasing data')
    check.one('phased.hap2', 'join phasing data')
    check.one('training.site', 'run annotate.gatk.training')
}


# Full GATK annotation pipeline. Writes out the table expected as
# input for SCAN2 call_mutations.
#
# IMPORTANT: mutsig annotations are NOT added here because they can
# be quite slow. These are better handled in parallel in chunks.
make.gatk.table <- function(mmq60, mmq1, phasing, single.cell, bulk, quiet=FALSE) {
    gatk <- read.gatk.table.2sample(mmq60, single.cell, bulk, region=NULL, quiet=quiet)
    annotate.gatk(gatk, add.mutsig=FALSE)
    annotate.gatk.lowmq(gatk, mmq1, single.cell, bulk, region=NULL, quiet=quiet)
    annotate.gatk.phasing(gatk, phasing)
    annotate.gatk.training(gatk, single.cell, bulk)
    check.gatk(gatk)
    gatk
}
