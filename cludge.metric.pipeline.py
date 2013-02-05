'''
This is going to be messy.

String together vcftools, plink, and python to produce a set of QC metrics on a
vcf.
'''

import subprocess as sp

devnull = open('/dev/null', 'w')
vcfFile = './data/ATLAS.merged.ontarget.chr22.SM.vcf'
sampleKeyFile = '/nas40t2/tsExoms/Sample Info/Scharf_master_sequencing_sample_key.csv'
pedFile = '/nas40t2/tsExoms/Sample Info/updated_Scharf_Pedigree.csv'
intermedBase = 'temp'


### get a MAF distribution
#mafCall = ['vcftools', '--vcf', vcfFile, '--freq', '--out', intermedBase] 
##print ' '.join(mafCall)
#sp.call(mafCall)
#
### get the ts/tv ratio
#tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', '1' , '--out', intermedBase]
##print ' '.join(tstvCall)
#sp.call(tstvCall)
#
### get the missingness rates
#missCall = ['vcftools', '--vcf', vcfFile, '--missing', '--out', intermedBase]
##print ' '.join(missCall)
#sp.call(missCall)
#
#### mendelian errors using plink
#
### convert vcf to plink formatted files
#convertCall = ['vcftools', '--vcf', vcfFile, '--plink', '--recode', '--out',
#    intermedBase]
#sp.call(convertCall)

## store mapping of SM to RG
rg_sm = dict()
fconn = open(sampleKeyFile)
header = fconn.next()
for line in fconn:
   line = line.split(',')
   rg_sm[line[1]]  = line[0]
fconn.close()

## add pedigree information to plink ped file file
mothers = dict()
fathers = dict()
families = dict()
sexes = dict()
with open(pedFile) as fconn:
    header = enumerate(fconn.next().split(','))
    index = dict( ( (x[1], x[0]) for x in header ) )
    for line in fconn:
        line = line.split(',')
        if rg_sm.get(line[index['Mother']]):
            mother = rg_sm[ line[index['Mother']] ]
        else:
            mother = line[index['Mother']]

        if rg_sm.get(line[index['Father']]):
            father = rg_sm[ line[index['Father']]]
        else:
            father = line[index['Father']]
        
        if mother=='Not Sequenced':
            mother = '0'
        if father=='Not Sequenced':
            father = '0'

        fathers[line[index['SM_tag (Local_ID)']]] = father
        mothers[line[index['SM_tag (Local_ID)']]] = mother

        if rg_sm.get(line[index['Individ']]):
            family = line[index['Family']]
        else:
            family = line[index['SM_tag (Local_ID)']] 

        families[line[index['SM_tag (Local_ID)']]] = family

        sexes[line[index['SM_tag (Local_ID)']]] = line[index['Sex']]


newPlinkFile = intermedBase + '.pedigree.ped'
newPlinkConn = open(newPlinkFile, 'w')
with open(intermedBase + '.ped') as plinkFile:
    for line in plinkFile:
        line = line.strip().split('\t')
        iid = line[1]
        if fathers.get(iid):
            line[2] = fathers[iid]
        if mothers.get(iid):
            line[3] = mothers[iid]
        if families.get(iid):
            line[0] = families[iid]
        if sexes.get(iid):
            line[4] = sexes[iid] 
        ## add a fake phenotype
        line[5] = '2'
        print >> newPlinkConn, '\t'.join(line)
newPlinkConn.close()

## calculate mendelian errors
mendelCall = ['plink',
    '--ped',
    newPlinkFile,
    '--map',
    intermedBase + '.map',
    '--mendel',
    '--allow-no-sex',
    '--out',
    intermedBase]
sp.call(mendelCall)#, stdout = devnull, stderr = devnull)



devnull.close()
