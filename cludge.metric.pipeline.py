'''
This is going to be messy.

String together vcftools, plink, and python to produce a set of QC metrics on a
vcf.
'''

import subprocess as sp


vcfFile = './data/ATLAS.merged.ontarget.chr22.SM.vcf'
intermedBase = 'temp'


## get a MAF distribution
mafCall = ['vcftools', '--vcf', vcfFile, '--freq', '--out', intermedBase] 
#print ' '.join(mafCall)
sp.call(mafCall)

## get the ts/tv ratio
tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', '1' , '--out', intermedBase]
#print ' '.join(tstvCall)
sp.call(tstvCall)

## get the missingness rates
missCall = ['vcftools', '--vcf', vcfFile, '--missing', '--out', intermedBase]
#print ' '.join(missCall)
sp.call(missCall)

## mendelian errors using plink



