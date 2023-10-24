# Number of haploid basepairs (in billions) per genome. The default
# value of 5.845001134 corresponds to AUTOSOMES as determined by GRCh37


get.gbp.by.genome <- function(object) {
    if (object@genome.string == 'hs37d5') {
        # 93 contigs includes unplaced; see http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics.
        total <- 3137161264
        chrx <- 155270560
        chry <- 59373566
        chrm <- 16571
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.845001134
    } else if (object@genome.string == 'hg38') {
        # 455 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 3209286105
        chrx <- 156040895 
	chry <- 57227415
        chrm <- 16569
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.992002452
    } else if (object@genome.string == 'mm10') {
        # 66 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 2730871774
        chrx <- 171031299
        chry <- 91744698
        chrm <- 16299
        return((total - chrx - chry - chrm)*2 / 1e9) # = 4.936158956
    } else {
        warn(paste('gbp not yet implemented for genome', object@genome.string))
        warn("the mutation burden for this analysis is a placeholder!")
        warn("DO NOT USE!")
        # hopefully returning a negative number will alert people that something
        # has gone wrong so they don't ignore the warning messages above
        return(-1)
    }
}


setGeneric("compute.mutburden", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE)
        standardGeneric("compute.mutburden"))
setMethod("compute.mutburden", "SCAN2", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE) {
    check.slots(object, c('call.mutations', 'depth.profile'))

    muttypes <- c('snv', 'indel')
    object@mutburden <- setNames(lapply(muttypes, function(mt) {
        # [2] is the maximum burden; the minimum burden [1] is almost always ~0
        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
        sfp <- object@static.filter.params[[mt]]
    
        # germline sites
        g <- object@gatk[resampled.training.site == TRUE & muttype == mt]

        # (single cell  x  bulk) depth table
        dptab <- object@depth.profile$dptab
        dptab <- dptab[1:min(max(g$dp)+1, nrow(dptab)),]

        # these computations rely on there being a reasonably large number
        # of germline sites tested. even 100 is very few; we expect more like
        # 100,000.
        if (nrow(g) < 100) {
            warning(paste('only', nrow(g), 'resampled germline', mt, 'sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites'))
            ret <- data.frame(
                ncalls=NA,
                callable.sens=NA,
                callable.bp=NA
            )[c(1,1,1),]
        } else {
            
			# trying to see these warnings
			summary(warnings())
			

			#attempting to load scan2 to try saving the object so I can actually look at it :(
			save(object, "test_scan2_pre_cut.rda")
			# somatic sites
            s <- object@gatk[pass == TRUE & muttype == mt]

            # Break data into 4 quantiles based on depth, use the middle 2 (i.e.,
            # middle 50%) to reduce noise caused by very low and very high depth.
            q=4
            qbreaks <- quantile(g$dp, prob=0:q/q)
    
            # s also uses g-based depth quantiles
			print("printing s")
			print(s)
			print("printing s$dp")
			print(s$dp)
            tryCatch(s$dpq <- cut(s$dp, qbreaks, include.lowest=T, labels=F,
					error = function(e) {
					message("cut to generate quantiles failed, is your rdata object generating a sensible dp column?")
			})
            s$dpq[s$dpq==3] <- 2 # merge 25-75% into a single bin
            g$dpq <- cut(g$dp, qbreaks, include.lowest=T, labels=F)
            g$dpq[g$dpq==3] <- 2

            # select the subset of the depth profile passing the bulk depth requirement
            # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
            rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
            rowqs[rowqs==3] <- 2
    
            qstouse <- c(1,2,4)
            s <- s[dpq %in% qstouse]
            g <- g[dpq %in% qstouse]
    
            # this data.frame has 1 row for each quantile. the second row (=middle 50%)
            # is ultimately what we're interested in, but having the other calculations
            # around can also be interesting.
            ret <- data.frame(
                ncalls=sapply(qstouse, function(q) sum(s[dpq == q]$pass, na.rm=TRUE)),
                callable.sens=sapply(qstouse, function(q) mean(g[bulk.dp >= sfp$min.bulk.dp & dpq == q]$resampled.training.pass, na.rm=TRUE)),
                callable.bp=sapply(split(dptab[,-(1:sfp$min.bulk.dp)], rowqs), sum)
            )
        }

        # "callable" means:
        # Sensitivity estimates only from germline training sites with the same
        # depth cutoffs as somatic candidates. Detailed depth tables will be used
        # to ensure extrapolation to the rest of the genome is equitable.
        ret$callable.burden <- ret$ncalls / ret$callable.sens
        # dividing by 2 makes it haploid gb
        ret$rate.per.gb <- ret$callable.burden / ret$callable.bp * 1e9/2
        ret$burden <- ret$rate.per.gb * gbp.per.genome
        ret$somatic.sens <- ret$ncalls / ret$burden
        ret$pre.genotyping.burden <- pre.geno.burden
        ret
    }), muttypes)

    object
})
