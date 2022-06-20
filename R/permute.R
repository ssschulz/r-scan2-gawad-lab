desired.perms <- as.integer(args[1])
sample <- args[2]
inrda <- args[3]
outfile <- args[4]
mutclass <- args[5]
callable.file <- args[6]
genome.file <- args[7]
seed.base <- as.integer(args[8])*1e6

if (!(mutclass %in% c('snv','indel')))
    stop('argument 5 must be either "snv" or "indel"')

if (file.exists(outfile))
    stop(paste('output file', outfile, 'already exists, please delete it first'))

if (!file.exists(callable.file))
    stop("BED file containing callable regions does not exist")

if (!file.exists(genome.file))
    stop("genome file containing genome coordinates does not exist")

cat('got seed.base and multiplied by 1e6', seed.base/1e6, '->', seed.base, '\n')


suppressMessages(library(bedtoolsr))
suppressMessages(library(scan2))
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.1000genomes.hs37d5))


# This function uses bedtools to randomly select positions (bedtools shuffle)
# within the region defined by 'callable'.
#
# n.sample - the number of permutations to create
# callable - a BED file with intervals representing valid regions to place
#          permuted mutations.
# genome - a BEDtools genome file, which lists one chromosome and its length
#          per line, separated by a tab.
bedtools.permute <- function(n.sample, genome, callable, seed) {
    g <- fread(genome)  # only read this to get a valid chromosome name

    real.n.sample <- n.sample
    n.sample <- n.sample*1.05 # add 5% to allow removal of positions < 50 bp
    if (n.sample > g[1,2][[1]])
        stop(paste('n.sample cannot exceed', g[1,2][[1]], 'due to current code limitations. '))

    # dummy data frame of single base positions to shuffle. the positions are ignored.
    tmpmuts <- cbind(rep(g[1], n.sample), 1:n.sample, 2:(n.sample+1))
    
    perms <- bt.shuffle(i=tmpmuts, g=genome, incl=callable, seed=seed, noOverlapping=TRUE)[,1:2]

    colnames(perms) <- c('chr', 'pos')
    # beds are 0-indexed, getSeq is 1-indexed. this matters in the case where
    # the bed returns position 0, which will cause getSeq to throw an error.
    # maybe affects some other corner cases as well.
    perms <- perms[perms$pos >= 50,]
    perms$chr <- as.character(perms$chr)
    perms <- head(perms, real.n.sample)
    if (nrow(perms) < n.sample)
        stop(paste('failed to generate enough sites. requested', n.sample, 'got', nrow(perms)))

    if (any(duplicated(paste(perms$chr, perms$pos))))
        stop('duplicate permutations found despite noOverlapping=TRUE')

    perms
}


################################################################################
# SNV PERMUTATION GENERATOR
################################################################################

# permute SNVs and preserve the mutation signature
# n.sample - this is the number of mutations to randomly permute around
#            the genome PRIOR to downsampling.
#            this number should be quite high, especially if the mutation
#            set has many mutated bases at CpGs.
#            XXX: eventually, it would be nice to use a smaller n.sample
#                 within a loop that could add more samples as necessary.
make.snv.perms.helper <- function(muts, genome, callable, seed, n.sample=5e4) {
    perms <- bedtools.permute(muts, n.sample=n.sample, genome=genome, callable=callable, seed=seed)
str(perms)

    # Get the reference base at each permutation position. Could save time here
    # and get the trinucleotide context, but oh well.
    perms$refnt <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5,
        names = perms$chr, start=perms$pos, end=perms$pos, as.character=TRUE)

    # rarely, a locus in the callable region can be an N
    # this can occur, e.g., if a stray single N is in the sequence or
    # because the ends of reads near an N gap can cover a few Ns despite
    # not matching them.
    # in any case, get rid of these.
    perms <- perms[perms$refnt != 'N',]

    # corner case:
    # issue arises when there are no input mutations starting from
    # a particular refnt. then there is no row in mutprobs for that
    # refnt and the subsequent subset will fail.  even if that subset
    # succeeded, there would be 0 probability for all altnts, and
    # sample will then fail.
    # if that is the case, just remove those permuted sites ahead of
    # time. this will only occur for samples with very small numbers
    # of mutations.
    perms <- perms[perms$refnt %in% unique(muts$refnt),]

    # get table of refnt > altnt probabilities to make sampling more efficient
    mat <- table(muts$refnt, muts$altnt)
    mutprobs <- mat/rowSums(mat)

    # Used to delete this, but it has been helpful to keep for when bugs pop up.
    print(table(perms$refnt))
    cat('mutprobs\n')
    print(mutprobs)

    perms$altnt <- apply(mutprobs[perms$refnt,,drop=FALSE], 1, function(row)
        sample(x=names(row), size=1, replace=FALSE, prob=row))
    perms <- get.3mer(perms)

    real.muts <- table(muts$type.and.ctx)
    perm.muts <- table(perms$type.and.ctx)[names(real.muts)]
    # how many permutation sets can we get from this sampling?
    limits <- floor(perm.muts/real.muts)
    print(cbind(real=real.muts, perm=perm.muts, ratio=limits))
    cat('top limiting factors\n')
    print(sort(limits))
    k <- min(limits)
    if (k < 1)
        stop("n.sample too low: unable to complete any permutations")
    
    # just randomly reorder (this isn't necessary, bt.shuffle is already unordered
    perms <- perms[sample(nrow(perms), size=nrow(perms), replace=FALSE),]
    # take the first N of each type and ctx. since order is random, this is equivalent
    # to selecting a random subset of each SBS channel.

    list(k=k, perms=do.call(rbind, lapply(names(real.muts), function(mt) {
        n.real <- real.muts[mt]
        head(perms[perms$type.and.ctx == mt,], k*n.real)
    })))
}




################################################################################
# INDEL PERMUTATION GENERATORS
################################################################################


# Updated version of classify.indels/muts to avoid the behavior of
# SigProfilerMatrixGenerator when there are duplicate indels.
classify.indels <- function(df, sample.name='dummy', save.plot=F, auto.delete=T, chrs=c(1:22,'X')) {
    ret <- classify.muts(df=df, spectype='ID', sample.name=sample.name, save.plot=save.plot,
        auto.delete=auto.delete, chrs=chrs)
    # remove transcribed strand information
    ret$muttype <- substr(ret$muttype, 3, 11)
    ret
}


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

    # SigProfilerMatrixGenerator doesn't classify duplicated mutations
    # from the same sample, it throws errors instead. It also will not
    # detect duplicates if they are not adjacent in the file.  If the
    # duplicate is not detected in this way, it causes the bug where
    # the final newdf dataframe does not match the input df.
    # To circumvent all these headaches: just remove duplicates up front.
    mutid <- paste(s$chr, s$pos, s$refnt, s$altnt)
    dupmut <- duplicated(mutid)
    cat("Removing", sum(dupmut), "/", nrow(s), "duplicated mutations before annotating\n")
    s <- s[!dupmut,]

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

    if (!all(df$chr == newdf$chr))
        stop('df and newdf do not perfectly correspond: df$chr != newdf$chr')
    if (!all(df$pos == newdf$pos))
        stop('df and newdf do not perfectly correspond: df$pos != newdf$pos')

    if (auto.delete)
        unlink(spmgd, recursive=TRUE)
    df$muttype <- newdf$iclass
    df
}


# For k=3, got ~100 of the rarest indel classes after removing
# sites in blacklist or within 200 bp of each other.
make.indel.perms.helper <- function(muts, spectrum, genome, callable, seed, k=1/10) {
    # The values below are sufficient to get about 100 of each indel type
    # however, getting the rare indel types is not as important for permutations
    # as it was for the original sensitivity testing. this is because we are
    # trying to match the indels called in the sample, in which there are very
    # few rare indel types.

    # this is an attempt to manually tune random indel generation so
    # that we can finish in reasonable time.
    # it will be better to do this automatically, as is done for SNVs
    ndels=10000000*k
    nins=500000*k
    nrins=5000000*k

    n.sample <- ndels+nins+nrins
    cat('generating', n.sample, 'candidates at a time\n')
    perms <- bedtools.permute(muts, n.sample=n.sample, genome=genome, callable=callable, seed=seed)
    cat('memory after generating permutations:\n')
    print(gc())


    # For DELETIONS
    cat('generating random deletions..\n')
    del.chrs <- perms$chr[1:ndels]
    del.locs <- perms$pos[1:ndels]
    del.lens <- 1+rnbinom(length(del.locs), mu=3, size=100)
    # XXX: AUTO del.lens <- sample(1:100, size=ndels, prob=del.ldist, replace=TRUE)
    del.refnts <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5,
        names=del.chrs, start=del.locs-1, end=del.locs+del.lens-1, as.character=TRUE)
    deld <- data.frame(chr=del.chrs, pos=del.locs-1, refnt=del.refnts, altnt=substr(del.refnts,1,1),
        stringsAsFactors=FALSE)
    cat('memory after generating deletions:\n')
    print(gc())
str(deld)


    # For INSERTIONS
    cat('generating random insertions..\n')
    bases <- c('A','C','G','T')
    ins.chrs <- perms$chr[ndels + (1:nins)]
    ins.locs <- perms$pos[ndels + (1:nins)]
    ins.lens <- 1+rnbinom(length(ins.locs), mu=3, size=100)
    # XXX: AUTO ins.lens <- sample(1:100, size=ndels, prob=ins.ldist+rins.ldist, replace=TRUE)
    ins.refnts <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5,
        names=ins.chrs, start=ins.locs, end=ins.locs, as.character=TRUE)
    ins.altnts <- sapply(ins.lens, function(l)
        paste0(sample(bases, l, replace=T), collapse=''))
    ind <- data.frame(chr=ins.chrs, pos=ins.locs, refnt=ins.refnts,
        altnt=paste0(ins.refnts, ins.altnts),
        stringsAsFactors=FALSE)
    cat('memory after generating insertions:\n')
    print(gc())
str(ind)

    # Repetitive INSERTIONS
    # no point in trying to generate repetitive elements perfectly.
    # just do something that generates them at much higher rates
    # than the above code; then downsample later.
    # heuristic:
    #  choose a position, unit size and #units
    #  extract 'unit' reference bases *after* the base at the chosen
    #  position, repeat those bases #units times
    cat('generating random insertions in repetitive sites..\n')
    rins.chrs <- perms$chr[ndels+nins + (1:nrins)]
    rins.locs <- perms$pos[ndels+nins + (1:nrins)]
    unit <- sample(5, size=length(rins.locs), replace=T)  # uniform on 1..5
    nunits <- sample(4, size=length(rins.locs), replace=T) # uniform on 1..4
    # get sequence for:
    #   [ pos | <-- unit --> ]
    #           ^^^^^^^^^^^^----- replicate this stretch nunits times
    # ..do this enough times that the replicated unit will intersect a
    # repetitive region of the same size.
    nts <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5,
        names=rins.chrs, start=rins.locs, end=rins.locs+unit, as.character=TRUE)
    repli <- substr(nts,2,nchar(nts))
    refnts <- substr(nts, 1, 1)
    extension <- sapply(1:length(rins.locs), function(i)
        paste0(rep(repli[i], nunits[i]), collapse=''))
    rind <- data.frame(chr=rins.chrs, pos=rins.locs, refnt=refnts,
        altnt=paste0(refnts, extension),
        stringsAsFactors=FALSE)
    cat('memory after generating repetitive insertions:\n')
    print(gc())
str(rind)


    cat('combining all random mutations..\n')
    d <- rbind(ind, deld, rind)
    # Don't allow any Ns in the reference or alternate sequence
    d <- d[!grepl('N', d$refnt) & !grepl('N', d$altnt),]
    d <- classify.indels(d)
cat(nrow(ind), '\n')
cat(nrow(deld), '\n')
cat(nrow(rind), '\n')
print(table(is.na(d$muttype)))
print(d[is.na(d$muttype),])
    cat('memory after classify.indels:\n')
    print(gc())


    # now that we have a table of random indels, figure out how
    # many permutation sets that match the provided spectrum can
    # be taken from the table.
    real.muts <- spectrum
    perm.muts <- table(d$muttype)[names(real.muts)]
    # how many permutation sets can we get from this sampling?
    limits <- floor(perm.muts/real.muts)
    print(cbind(real=real.muts, perm=perm.muts, ratio=limits))
    cat('top limiting factors\n')
    print(sort(limits))
    k <- min(limits)
    if (k < 1)
        stop("k is too low: unable to complete any permutations")
    
    # just randomly reorder (this isn't necessary, bt.shuffle is already unordered
    d <- d[sample(nrow(d), size=nrow(d), replace=FALSE),]
    cat('memory after resampling d:\n')
    print(gc())

    # take the first N of each type and ctx. since order is random, this is equivalent
    # to selecting a random subset of each SBS channel.
    # these will be chopped up and turned into valid permutation sets by the caller.
    return(
        list(k=k, perms=do.call(rbind, lapply(names(real.muts), function(mt) {
            n.real <- real.muts[mt]
            head(d[d$muttype == mt,], k*n.real)
        })))
    )
}
    


# basic idea: each call to make.perms.helper creates a large random
# sample of SNVs, typically enough to create several permutation sets.
# just continue to call the helper function until all desired perms
# are created.
# this allows us to bypass the small, constant overhead of bt.shuffle,
# which is 5-10 seconds even when the BED file to be shuffled has as
# few as 10 lines.
# example:
#   using n.samples=5 million, 218 permutations of 5657-Oligo-7 (1023 sSNVs)
#   can be solved in a single helper() call. (4 minute runtime)
#   using n.samples=10 million, 6168 permutations of 1278-Oligo-5 (56 sSNVs)
#   can be solved per call. (8 minute runtime)
#
# max RAM usage is ~4GB for n.samples=10 million or k=1/5
make.perms <- function(muts, genome, callable, mutclass=c('snv','indel'), desired.perms=10000, ...) {
    mutclass <- match.arg(mutclass)
    permuted.muts <- NULL
    i <- 1
    total.solved <- 0
    seeds.used <- c()

    # drop out early if there are no mutations
    if (nrow(muts) == 0) {
        return(list(total.solved=desired.perms,
            seeds.used=seed.base+1,                         # whatever
            muts=muts, raw.perms=NULL,
            perms=lapply(1:desired.perms, function(i) NULL) # list of NULLs
        ))
    }

    # SNV generator does this internally, though has an ignored 'spectrum' argument
    # for compatible calling.
    if (mutclass == 'indel') {
        cat('Making INDEL permutations\n')
        spectrum.to.duplicate <- plot.indel(iclass=muts$muttype, make.plot=FALSE, reduce.to.id83=FALSE)
    } else {
        cat('Making SNV permutations\n')
    }

    while (desired.perms - total.solved > 0) {
        this.seed <- seed.base+i
        # list of seeds that, for some reason, cause a floating point exception
        # in bedtools shuffle. this only happened for one seed at depths 25 and 30.
        if (this.seed %in% c(329000001))
            this.seed <- this.seed+1
        cat('iteration', i, 'remaining to solve', desired.perms - total.solved, 'seed', this.seed, '\n')
        if (mutclass == 'indel') {
            print(system.time(
                ret <- make.indel.perms.helper(muts=muts, spectrum=spectrum.to.duplicate,
                    genome=genome, callable=callable, seed=this.seed, ...)
            ))
        } else {
            print(system.time(
                ret <- make.snv.perms.helper(muts=muts,
                    genome=genome, callable=callable, seed=this.seed, ...)
            ))
        }
        cat('memory after make.XXX.perms.helper:\n')
        print(gc())
        i <- i+1
        seeds.used <- c(seeds.used, this.seed)
        total.solved <- total.solved + ret$k
        permuted.muts <- rbind(permuted.muts, ret$perms)
        # overshooting the requested perms is fine
    }

    # chop up the mutations into individual permutation sets
    # sadly, SNVs and indels save signature components in different columns
    if (mutclass == 'snv') {
        muttype.counts <- table(muts$type.and.ctx)
        muts.by.muttype <- split(permuted.muts, permuted.muts$type.and.ctx)
    } else if (mutclass == 'indel') {
        muttype.counts <- table(muts$muttype)
        muts.by.muttype <- split(permuted.muts, permuted.muts$muttype)
    }
    all.muttypes <- names(muttype.counts)
    str(muts.by.muttype)
    print(all.muttypes)
    print(muttype.counts)
    perms <- lapply(1:desired.perms, function(i) {
        do.call(rbind, lapply(all.muttypes, function(mt)
            muts.by.muttype[[mt]][(1 + (i-1)*muttype.counts[mt]):(i*muttype.counts[mt]),]))
    })
    list(total.solved=total.solved, seeds.used=seeds.used, muts=muts, raw.perms=permuted.muts, perms=perms)
}





muts <- get(load(inrda))
muts <- muts[muts$sample == sample,]
cat("Found", nrow(muts), "mutations of type", mutclass, "for sample", sample, '\n')

if (mutclass == 'snv') {
    system.time(permdata <- make.perms(
        muts=muts,
        mutclass=mutclass,
        genome=genome.file,
        callable=callable.file,
        desired.perms=desired.perms, n.sample=1e6))
        # real usage
        #n.sample=1e6))
        # smaller for testing
        #n.sample=1e4))
} else { 
    system.time(permdata <- make.perms(
        muts=muts,
        mutclass=mutclass,
        genome=genome.file,
        callable=callable.file,
        desired.perms=desired.perms, k=1/5))
        # real usage
        # N.B.: potential duplicate positions are removed on each loop, the number of
        # which depends on k. higher values of k mean fewer loops and thus lower likelihood
        # of duplicate positions, which may be undesirable.
        # duplicate rates are generally low and it is debatable whether they should be
        # removed at all--they tend to represent regions with local sequence contexts
        # amenable to creating certain types of indels and thus may also be where
        # mutations of the same type are also more frequent. Note there is no removal of
        # duplicate permuted indels across samples.
        #k=1/5))
        # smaller for testing
        #k=1/100))
}

cat('final memory usage:\n')
print(gc())

save(sample, inrda, permdata, genome.file, callable.file, seed.base, file=outfile, compress=FALSE)
