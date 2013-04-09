## there is an array of input files:
##  * vcftools imiss
##  * vcftools lmiss
##  * vcftools frq
##  * vcftools TsTv
##  * plink imendel
##  * plink lmendel

library(ggplot2)

atlas.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_data/atlas_exome_chrid_sm.vcf'
freebayes.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_data/freebayes_03-25_exome_minQ.vcf'
gatk.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_data/gatk_05-05_exome.vcf'
consensus.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_data/ts.exome.consensus.vcf'

bases <- c(atlas.base, freebayes.base, gatk.base, consensus.base)

## locus missingness
strip.lmiss <- function(fname) {
    dat <- read.table(fname, header=T)
    return( dat$F_MISS )
}

strip.branch.name <- function(fname) {
    tmp <- tail(unlist(strsplit(fname,'/')), 1)
    name <- head(unlist(strsplit(tmp, '\\.')), 1)
    return(name)
}

##missingness spectra
lmiss <- c()
for (i in seq_along(bases)) {
    branch.name <- strip.branch.name(bases[i])
    file.path <- tmpFile <- paste(bases[i], '.lmiss', sep='')
    lmiss.dat <- read.table(file.path, header=T)
    lmiss <- c(lmiss, list(lmiss.dat$F_MISS))
}
miss.df <- rbind(data.frame(dataset='ATLAS', missf=1.0 - unlist(lmiss[1])),
        data.frame(dataset='Freebayes', missf=unlist(lmiss[2])),
        data.frame(dataset='GATK', missf=unlist(lmiss[3])),
        data.frame(dataset='Consensus', missf=unlist(lmiss[4]))
        )
plt <- ggplot(miss.df, aes(x=missf, fill=dataset)) +
  geom_histogram(colour="black", position="dodge") +
  scale_fill_manual(values=c("blue","green","red","orange"))
show(plt)

mafs <- c()
## locus MAF spectrum
for (i in seq_along(bases)) {
    branch.name <- strip.branch.name(bases[i])
    file.path <- tmpFile <- paste(bases[i], '.frq', sep='')
    maf.dat <- read.table(file.path, header=T)
    mafs <- c(mafs, list(maf.dat$MAF))
}
maf.df <- rbind(data.frame(dataset='ATLAS', maf=0.5-unlist(mafs[1])),
        data.frame(dataset='Freebayes', maf=unlist(mafs[2])),
        data.frame(dataset='GATK', maf=unlist(mafs[3])),
        data.frame(dataset='Consensus', maf=unlist(mafs[4]))
        )
DF$dataset <- as.factor(maf.df$dataset)
plt <- ggplot(DF, aes(x=maf, fill=dataset)) +
  geom_histogram(colour="black", position="dodge") +
  scale_fill_manual(values=c("blue","green","red","orange"))
show(plt)

## gstenome Ts/Tv ratio
atlas.tstv <- 68980.0 / 22988.0
freebayes.tstv <- 37937.0 / 31580.0
gatk.tstv <- 68243.0 / 27671.0
consensus.tstv <- 16497.0 / 5131.0
tstv.dat <- data.frame(c(atlas.tstv, freebayes.tstv, gatk.tstv,
    consensus.tstv))

x11()
barplot(tstv.dat[,1], names.arg=c('atlas', 'freebayes', 'gatk', 'consensus'),
    main = 'Ts/Tv')

## number of variants
atlas.m <- 91968
freebayes.m <- 69517
gatk.m <- 95914
consensus.m <- 21628
n.dat <- data.frame(c(atlas.m,freebayes.m,gatk.m,consensus.m))
x11()
barplot(n.dat[,1], names.arg=c('atlas', 'freebayes', 'gatk', 'consensus'),
    main='Number of Variants')


## locus mendelian inconsistencies
lmendel.n <- c()
names <- c()
for (i in seq_along(bases)) {
    branch.name <- strip.branch.name(bases[i])
    names <- c(names, branch.name)
    lmendelFile <- paste(bases[i], '.lmendel', sep='')
    lmendelDat <- read.table(lmendelFile, header=T)
    lmendel.n <- c(lmendel.n, sum(as.numeric(lmendelDat$N)))
}
x11()
barplot(data.frame(lmendel.n)[,1], names.arg=c('atlas', 'freebayes', 'gatk',
'consensus'),
    main='Mendelian Inconsistencies' )

