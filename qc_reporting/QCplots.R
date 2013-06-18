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


## read in data for loose consensus base
qual.plot <- function(fpath, name){
    dat <- read.table(fpath, na.strings='-', header=T,
      colClasses=c(rep('character',4), rep('numeric',3)))
    plt <- ggplot(dat, aes(x=log10(FQ), y=log10(GQ))) + geom_point(alpha=0.4) +
      geom_abline(aes(intercept=0, slope=1)) +
      labs(title=paste(name, "QUAL score comparison"), x="log freebayes QUAL",
          y="log GATK qual")
    show(plt)
}


### INPUTS ###
##############

atlas.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/atlas_exonBed_v1.4.3_allsamples_flat/atlas_exome_bed_v1.4.3_allsamples_flat' 
freebayes.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/freebayes_exonBed_filter/freebayes_03-25_exome_bed'
gatk.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/gatk_exonBed_filter/gatk_05-05_exome_bed'
loose.consensus.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/loose_exonBed_v1.4.3_flat/consensus.loose.exonBed.atlas1.4.3.allsamples.flat'
strict.consensus.base <-
'/nas40t0/vasya/consensus_call/galaxy.consensus/ts_exomes_qc/strict_exonBed_v1.4.3_flat/consensus.strict.exonBed.atlas1.4.3.allsamples.flat'

bases <- c(atlas.base, freebayes.base, gatk.base, strict.consensus.base,
  loose.consensus.base)

base.names <- c('atlas', 'freebayes', 'gatk', 'strict', 'loose')


## calculate missingness spectra
lmiss <- c()
total.loci <- 0
for (i in seq_along(bases)) {
    file.path <- tmpFile <- paste(bases[i], '.lmiss', sep='')
    lmiss.dat <- read.table(file.path, header=T)
    lmiss <- c(lmiss, list((1.0 - lmiss.dat$F_MISS)))
    total.loci <- total.loci + length(lmiss.dat$F_MISS)
}
miss.df <- rbind(data.frame(dataset=base.names[1], missf=unlist(lmiss[1])),
        data.frame(dataset=base.names[2], missf=unlist(lmiss[2])),
        data.frame(dataset=base.names[3], missf=unlist(lmiss[3])),
        data.frame(dataset=base.names[4], missf=unlist(lmiss[4])),
        data.frame(dataset=base.names[5], missf=unlist(lmiss[5])))
 


## locus MAF spectrum
mafs <- c()
for (i in seq_along(bases)) {
    file.path <- paste(bases[i], '.frq', sep='')
    maf.dat <- read.table(file.path, header=T)
    mafs <- c(mafs, list(maf.dat$MAF))
}
maf.df <- rbind(data.frame(dataset=base.names[1], maf=unlist(mafs[1])),
        data.frame(dataset=base.names[2], maf=unlist(mafs[2])),
        data.frame(dataset=base.names[3], maf=unlist(mafs[3])),
        data.frame(dataset=base.names[4], maf=unlist(mafs[4])),
        data.frame(dataset=base.names[5], maf=unlist(mafs[5])))

## genome Ts/Tv ratio
tstv.dat <- vector(mode="list", length(base.names)) 
names(tstv.dat) <- base.names
for (i in seq_along(bases)) {
    file.path <- paste(bases[i], '.tstv.TsTv.summary', sep='')
    dat <- read.table(file.path, header=T)
    tstv.dat[base.names[i]] <- dat$COUNT[7] / dat$COUNT[8]
}



## locus mendelian inconsistencies and number of variants
lmendel.dat <- vector(mode="list", length(base.names))
names(lmendel.dat) <- base.names
mvar.dat <- vector(mode="list", length(base.names))
names(mvar.dat) <- base.names
for (i in seq_along(bases)) {
    file.path <- paste(bases[i], '.lmendel', sep='')
    dat <- read.table(file.path, header=T)
    lmendel.dat[base.names[i]] <- sum(as.numeric(dat$N))
    mvar.dat[base.names[i]] <- length(dat$SNP)
}




### PLOT GENERATION
## note we step through with the enter

## missingness
#plt2 <- ggplot(miss.df, aes(x=missf, color=dataset)) + stat_ecdf(geom="smooth")
#show(plt2)

#readline()
plt2 <- ggplot(miss.df, aes(x=missf, fill=dataset)) +
geom_histogram(binwidth=0.01) +
    facet_grid(dataset~.) + labs(title="Missingness Spectra",
        xlab="rate of missing genotypes")
show(plt2)

## regular histogram of missingness
#plt <- ggplot(miss.df, aes(x=missf, fill=dataset)) +
#  geom_histogram(colour="black", position="dodge",
#    aes(y=..count../ total.loci )) +
#  scale_fill_manual(values=c("blue","green","red","orange", "purple"))
#show(plt)

## MAF
#readline()
#plt <- ggplot(maf.df, aes(x=maf, fill=dataset)) +
#  geom_histogram(position="dodge", binwidth = 0.05) 
#show(plt)

readline()
plt <- ggplot(maf.df, aes(x=maf, fill=dataset)) +
geom_histogram(binwidth=0.005) +
    facet_grid(dataset~.) + labs(title="Minor Allele Frequency Distributions",
        x="maf")
show(plt)

### separate out individual data sets
#readline()
### this was an attempt to reweight the densities by number of variants
##plt <- ggplot(maf.df, aes(x=maf, color=dataset)) +
##    geom_density(data=subset(maf.df, dataset=='freebayes'), adjust=8,
##      weights=1/length(which(maf.df$dataset=='freebayes'))) +
##    geom_density(data=subset(maf.df, dataset=='atlas'), adjust=0.5,
##      weights=1/length(which(maf.df$dataset=='atlas'))) +
##    geom_density(data=subset(maf.df, dataset=='gatk'), adjust=0.5,
##      weights=1/length(which(maf.df$dataset=='gatk'))) +
##    geom_density(data=subset(maf.df, dataset=='strict'), adjust=1,
##      weights=1/length(which(maf.df$dataset=='strict'))) +
##    geom_density(data=subset(maf.df, dataset=='loose'), adjust=0.5,
##      weights=1/length(which(maf.df$dataset=='loose'))) 
#plt <- ggplot(maf.df, aes(x=maf, color=dataset)) +
#    geom_density(data=subset(maf.df, dataset=='freebayes'), adjust=8,
#binwidth=0.01) +
#    geom_density(data=subset(maf.df, dataset=='atlas'), adjust=0.5,
#binwidth=0.01) +
#    geom_density(data=subset(maf.df, dataset=='gatk'), adjust=0.5,
#binwidth=0.01) +
#    geom_density(data=subset(maf.df, dataset=='strict'), adjust=1,
#binwidth=0.01) +
#    geom_density(data=subset(maf.df, dataset=='loose'), adjust=0.5,
#binwidth=0.01) 
#show(plt)
#
### experimental ecdf representation
#readline()
#plt2 <- ggplot(maf.df, aes(x=maf, color=dataset)) + stat_ecdf(geom="smooth") +
#    labs(title="Minor Allele Frequency Cumulative Distribution Function")
#show(plt2)

## ts/tv
readline()
tstv.df <- data.frame(dataset=base.names, tstv=unlist(tstv.dat) )
plt <- ggplot(tstv.df, aes(y=tstv, x=dataset, fill=dataset)) + geom_bar() +
    labs(title="Ts/Tv Ratio", y="ts/tv ratio")
show(plt)

## number of variants
readline()
mvar.df <- data.frame(dataset=base.names, mvar=unlist(mvar.dat))
plt <- ggplot(mvar.df, aes(x=dataset, y=mvar, fill=dataset)) + geom_bar() +
    labs(title="Number of Variants", y="number of variants")
show(plt)

## mendel errors
readline()
lmendel.df <- data.frame(dataset=base.names,
lmendel.rate=unlist(lmendel.dat)/unlist(mvar.dat), lmendel=unlist(lmendel.dat))
plt <- ggplot(lmendel.df, aes(x=dataset, y=lmendel.rate, fill=dataset)) + geom_bar() +
    labs(title="Mendelian Error Rate", y="mendel error rate")
show(plt)


## QUAL distributions
readline()
qual.plot(paste(strict.consensus.base, '.INFO', sep=''), 'strict')

readline()
qual.plot(paste(loose.consensus.base, '.INFO', sep=''), 'loose')

## HWE
plot.new()
text(0.5, 0.5, 'HWE data not yet available')




