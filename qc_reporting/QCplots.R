## there is an array of input files:
##  * vcftools imiss
##  * vcftools lmiss
##  * vcftools frq
##  * vcftools TsTv
##  * plink imendel
##  * plink lmendel

### LIBRARIES ###
#################

library(ggplot2)


### FUNCTIONS ###
#################

## return vecotr of locus missingness
strip.lmiss <- function(fname) {
    dat <- read.table(fname, header=T)
    return( dat$F_MISS )
}

## return branch name
strip.branch.name <- function(fname) {
    tmp <- tail(unlist(strsplit(fname,'/')), 1)
    name <- head(unlist(strsplit(tmp, '\\.')), 1)
    return(name)
}


### INPUTS ###
##############

atlas.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/ts.exome.atlas' 
freebayes.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/ts.exome.freebayes'
gatk.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/ts.exome.gatk'
loose.consensus.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/ts.exome.consensus.loose'
strict.consensus.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/ts.exome.consensus.strict'

bases <- c(atlas.base, freebayes.base, gatk.base, strict.consensus.base,
loose.consensus.base)

base.names <- c('atlas', 'freebayes', 'gatk', 'loose', 'strict')


### calculate missingness spectra
#x11()
#lmiss <- c()
#total.loci <- 0
#for (i in seq_along(bases)) {
#    branch.name <- strip.branch.name(bases[i])
#    file.path <- tmpFile <- paste(bases[i], '.lmiss', sep='')
#    lmiss.dat <- read.table(file.path, header=T)
#    lmiss <- c(lmiss, list((1.0 - lmiss.dat$F_MISS)))
#    total.loci <- total.loci + length(lmiss.dat$F_MISS)
#}
#miss.df <- rbind(data.frame(dataset=base.names[1], missf=unlist(lmiss[1])),
#        data.frame(dataset=base.names[2], missf=unlist(lmiss[2])),
#        data.frame(dataset=base.names[3], missf=unlist(lmiss[3])),
#        data.frame(dataset=base.names[4], missf=unlist(lmiss[4])),
#        data.frame(dataset=base.names[5], missf=unlist(lmiss[5])))
# 
##plt <- ggplot(miss.df, aes(x=missf, fill=dataset)) +
##  geom_histogram(colour="black", position="dodge",
##    aes(y=..count../ total.loci )) +
##  scale_fill_manual(values=c("blue","green","red","orange", "purple"))
##show(plt)
#
#plt2 <- ggplot(miss.df, aes(x=missf, color=dataset)) + stat_ecdf(geom="smooth")
#show(plt2)
#
### locus MAF spectrum
#x11()
#mafs <- c()
#for (i in seq_along(bases)) {
#file.path <- tmpFile <- paste(bases[i], '.frq', sep='')
#    maf.dat <- read.table(file.path, header=T)
#    mafs <- c(mafs, list(maf.dat$MAF))
#}
#maf.df <- rbind(data.frame(dataset=base.names[1], maf=unlist(mafs[1])),
#        data.frame(dataset=base.names[2], maf=unlist(mafs[2])),
#        data.frame(dataset=base.names[3], maf=unlist(mafs[3])),
#        data.frame(dataset=base.names[4], maf=unlist(mafs[4])),
#        data.frame(dataset=base.names[5], maf=unlist(mafs[5])))
barplot(data.frame(lmendel.n)[,1], names.arg=base.names,
    main='Mendelian Inconsistencies' )
#
##DF$dataset <- as.factor(maf.df$dataset)
#plt <- ggplot(maf.df, aes(x=maf, fill=dataset)) +
#  geom_histogram(position="stack", binwidth = 0.1) 
#show(plt)
#
#plt2 <- ggplot(maf.df, aes(x=maf, color=dataset)) + stat_ecdf(geom="smooth")
#show(plt2)
#
#
#### gstenome Ts/Tv ratio
##atlas.tstv <- 68980.0 / 22988.0
##freebayes.tstv <- 37937.0 / 31580.0
##gatk.tstv <- 68243.0 / 27671.0
##consensus.tstv <- 16497.0 / 5131.0
##tstv.dat <- data.frame(c(atlas.tstv, freebayes.tstv, gatk.tstv,
##    consensus.tstv))
#
#tstv.df <- rbind(data.frame(dataset='ATLAS', tstv=0,
#        data.frame(dataset='Freebayes', tstv=0,
#        data.frame(dataset='GATK', tstv=0,
#        data.frame(dataset='Consensus', tstv=0
#)
#for (i in seq_along(bases)) {
#  branch.name <- strip.branch.name(bases[i])
#    file.path <- tmpFile <- paste(bases[i], '.frq', sep='')
#    maf.dat <- read.table(file.path, header=T)
#    
#
#}
#
##x11()pretty
##barplot(tstv.dat[,1], names.areg=c('atlas', 'freebayes', 'gatk', 'consensus'),
##    main = 'Ts/Tv')
##
#### number of variants
##atlas.m <- 10
##freebayes.m <- 10
##gatk.m <- 10
##consensus.m <- 10
##n.dat <- data.frame(c(atlas.m,freebayes.m,gatk.m,consensus.m))
##x11()
##barplot(n.dat[,1], names.arg=c('atlas', 'freebayes', 'gatk', 'consensus'),
##    main='Number of Variants'rplot(data.frame(lmendel.n)[,1],
#names.arg=base.names,
#    main='Mendelian Inconsistencies' )
#)
#
#
>>>>>>> qc_report
## locus mendelian inconsistencies
lmendel.n <- c()
names <- c()
for (i in seq_along(bases)) {
    lmendelFile <- paste(bases[i], '.lmendel', sep='')
    lmendelDat <- read.table(lmendelFile, header=T)
    lmendel.n <- c(lmendel.n, sum(as.numeric(lmendelDat$N)))
}
x11()
barplot(data.frame(lmendel.n)[,1], names.arg=base.names,
    main='Mendelian Inconsistencies' )

