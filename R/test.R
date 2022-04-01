testpipe <- function(use.big=FALSE) {
    if (!use.big) {
        # test parameters and files
        sc.sample.test <- 'h25'
        bulk.sample.test <- 'hunamp'
        mmq60.test <- system.file("extdata", "mmq60.tab.bgz", package="scan2")
        mmq1.test <- system.file('extdata', 'mmq1.tab.bgz', package='scan2')
        hsnps.test <- system.file('extdata', 'hsnps.tab.bgz', package='scan2')
        abfits.test <- system.file('extdata', 'fits.rda', package='scan2')
        sccigars.test <- system.file('extdata', 'h25_somatic_cigars.tab.bgz', package='scan2')
        bulkcigars.test <- system.file('extdata', 'hunamp_somatic_cigars.tab.bgz', package='scan2')
    } else {
        # much larger full chromosome and multi-sample dataset
        sc.sample.test <- '5823PFC-A'
        bulk.sample.test <- '5823-tempmusc-1b1_20170221-WGS'
        mmq60.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/5823_mmq60.tab.gz'
        mmq1.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/5823_mmq1.tab.gz'
        hsnps.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/5823_hsnps.tab.gz'
        abfits.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/5823PFC-A_fits.rda'
        sccigars.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/5823PFC-A_somatic_cigars.tab.gz'
        bulkcigars.test <- '/Users/jl/pta/pta_aging_reanalysis/scan2_dev/5823_data/dev/bulk_somatic_cigars.tab.gz'
    }




    #grs <- list(`1MB region 30M-31M`=GRanges(seqnames=22, ranges=IRanges(start=30000000, end=30999999)),
        #`full chr22`=GRanges(seqnames=22, ranges=IRanges(start=16000000, end=49999999)))
        #`full chr1`=GRanges(seqnames=1, ranges=IRanges(start=1, end=249999999)))

    #grs <- list( `full chr1`=GRanges(seqnames=1, ranges=IRanges(start=10e6, end=20e6)))


    grs <- list(`1MB region 30M-31M`=GRanges(seqnames=22, ranges=IRanges(start=30e6, end=30999999)),
                #`1MB region 31M-32M`=GRanges(seqnames=22, ranges=IRanges(start=31e6, end=31999999)),
                `1MB region 32M-33M`=GRanges(seqnames=22, ranges=IRanges(start=32e6, end=32999999)))
    lapply(1:length(grs), function(i) {
        cat('testing', names(grs)[i], '\n')
        gr <- grs[[i]]
        gt <- make.scan(single.cell=sc.sample.test, bulk=bulk.sample.test, genome='hs37d5', region=gr)
        cat('Before read.gatk:\n') ; print(gc(reset=TRUE))
        print(system.time(x1 <- read.gatk(gt, path=mmq60.test)))
        #print(x1)
        cat('After read.gatk:\n'); print(gc(reset=TRUE))
        print(system.time(y1 <- read.gatk.lowmq(x1, path=mmq1.test)))
        #print(y1)
        cat('After read.gatk.lowmq:\n') ; print(gc(reset=TRUE))
        print(system.time(z1 <- add.training.data(y1, path=hsnps.test)))
        #print(z1)
        cat('After add.training.data:\n') ; print(gc(reset=TRUE))
        print(system.time(zz1 <- resample.training.data(z1)))
        #print(zz1)
        cat('After resample.training.data:\n') ; print(gc(reset=TRUE))
        print(system.time(w1 <- add.ab.fits(zz1, path=abfits.test)))
        #print(w1)
        cat('After add.ab.fits:\n') ; print(gc(reset=TRUE))
        print(system.time(v1 <- compute.ab.estimates(w1)))
        #print(v1)
        cat('After compute.ab.estimates:\n') ; print(gc(reset=TRUE))
        print(system.time(r1 <- add.cigar.data(v1, sccigars.test, bulkcigars.test)))
        #print(r1)
        cat('After add.cigar.data:\n') ; print(gc(reset=TRUE))
        print(system.time(s1 <- compute.models(r1)))
        print(s1)
        cat('After compute.models:\n') ; print(gc(reset=TRUE))
        print(system.time(t1 <- compute.excess.cigar.scores(s1)))
        print(t1)
        cat('After compute.excess.cigar.scores:\n') ; print(gc(reset=TRUE))
        t1
    })
}
