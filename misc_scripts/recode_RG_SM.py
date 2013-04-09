'''
recode the plink tfam file, taking RG-SM IDs and mapping them to RG IDs
'''

import gzip
import os

BAMdir = '/nas40t2/tsExoms/BAM_BAI/ReadGroup/'
sm_to_rg = {}
for path, dir, files in os.walk(BAMdir):
    for file in files:
       if file.split('.')[-1]=='bam':
            for line in gzip.open(path+file):
                if line[0:3]=='@RG':
                    ids = line
                    break
                
            ids = ids.split('\t')
            for id in ids:
                if id[0:2]=='ID':
                    rg = id.split(':')[-1]
                elif id[0:2]=='SM':
                    sm = id.split(':')[-1]
            sm_to_rg[sm] = rg
            sm_to_rg[rg] = sm

oldvcfFile = open('ts_exomes_data/atlas_exome_chrid.vcf')
newvcfFile = open('ts_exomes_data/atlas_exome_chrid_sm.vcf', 'w')

line = oldvcfFile.next()
newvcfFile.write(line)
while not line.startswith('#CHROM'):
  newvcfFile.write(line)
  line = oldvcfFile.next()

newheader = []
for f in line.split():
  f = f.rstrip('.bam.snp.vcf')
  if sm_to_rg.get(f):
    newheader.append(sm_to_rg[f])
  else:
    newheader.append(f)
print >> newvcfFile, '\t'.join(newheader)


for line in oldvcfFile:
  newvcfFile.write(line)

