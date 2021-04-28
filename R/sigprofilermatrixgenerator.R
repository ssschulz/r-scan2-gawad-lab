# use sigprofilermatrixgenerator to classify mutations
# note that sigprofilermatrixgenerator always generates the most granular classification
# (e.g., SBS6144 for SNVs and ID415 for indels)
# mostly useful for indel classification, but also good for doing TSB on SNVs
# IMPORTANT BUG: the mutation classes are not always in the right order!!!
# this happens to not affect the SCAN* pipelines because the SCAN mutation
# order matches this order by chance. but please fix this later.
classify.muts <- function(df, spectype='SNV',
    sample.name='dummy', save.plot=F, auto.delete=T, chrs=1:22)
{
    if (nrow(df) == 0)
        return(df)

    recognized.spectypes <- c('SNV', 'ID')
    if (!(spectype %in% recognized.spectypes))
        stop(sprintf("unrecognized spectype '%s', currently only supporting %s",
            spectype, paste('"', recognized.spectypes, '"', collapse=', ')))

    require(SigProfilerMatrixGeneratorR)
    td <- paste0(tempdir())
    spmgd <- paste0(td, "/spmgr/")
    if (file.exists(spmgd) & auto.delete)
        unlink(spmgd, recursive=TRUE)
    dir.create(spmgd, recursive=TRUE)

    # convert this to use scansnv.df.to.vcf at some point
    # Write out the VCF
    ### From scansnv.to.vcf
    out.file <- paste0(spmgd, sample.name, '.vcf')
    f <- file(out.file, "w")
    vcf.header <- c("##fileformat=VCFv4.0", "##source=scansnv", 
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", 
            #sprintf("##contig=<ID=%s,length=%d>", fai[, 1], fai[,2]),
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", 
            "QUAL", "FILTER", "INFO", "FORMAT", sample.name), collapse = "\t"))
    writeLines(vcf.header, con = f)
    s <- df[!is.na(df$chr),]
    s <- do.call(rbind, lapply(chrs, function(chr) {
        ss <- s[s$chr == chr, ]
        ss[order(ss$pos), ]
    }))
    # will write, eg., position=7000000 as 7e6, which will
    # confuse sigprofilermatrixgenerator
    old.opt <- options('scipen')$scipen
    options(scipen=10000)
    writeLines(paste(s$chr, s$pos, s$dbsnp, s$refnt, s$altnt, 
        ".", "PASS", ".", "GT", "0/1", sep = "\t"), con = f)
    close(f)
    options(scipen=old.opt)

    mat <- SigProfilerMatrixGeneratorR(sample.name, 'GRCh37', spmgd, seqInfo=TRUE, plot=TRUE)

    # Read in the types
    annot.files <- paste0(spmgd, '/output/vcf_files/', spectype, '/', c(1:22,'X','Y'), "_seqinfo.txt")
    if (spectype == 'ID') {
        colclasses <- c(V2='character', V5='character', V6='character')
    } else if (spectype == 'SNV') {
        colclasses <- c(V2='character')
    }
    annots <- do.call(rbind, lapply(annot.files, function(f) {
        tryCatch(x <- read.table(f, header=F, stringsAsFactors=FALSE,
                colClasses=colclasses),
            error=function(e) NULL)
    }))
    if (spectype == 'ID') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'refnt', 'altnt', 'unknown')
        newdf <- plyr::join(df, annots[2:6], by = colnames(annots)[-c(1, 4, 7)])
    } else if (spectype == 'SNV') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'unknown')
        newdf <- plyr::join(df, annots[2:4], by = colnames(annots)[-c(1, 4, 5)])
    }

    if (save.plot) {
        plotfiles <- list.files(paste0(spmgd, '/output/plots/'), full.names=T)
        file.copy(plotfiles, '.')
    }

    if (auto.delete)
        unlink(spmgd, recursive=TRUE)
    if (!all(df$chr == newdf$chr))
        stop('df and newdf do not perfectly correspond: df$chr != newdf$chr')
    if (!all(df$pos == newdf$pos))
        stop('df and newdf do not perfectly correspond: df$pos != newdf$pos')
    df$muttype <- newdf$iclass
    df
}


# Convenience wrapper for older code
classify.indels <- function(df, sample.name='dummy', save.plot=F, auto.delete=T, chrs=1:22) {
    classify.muts(df=df, spectype='ID', sample.name=sample.name, save.plot=save.plot,
        auto.delete=auto.delete, chrs=chrs)
}
