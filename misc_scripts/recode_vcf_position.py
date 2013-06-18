'''
recode the chromosomal identifier in a vcf
'''

oldvcf = open('ts_exomes_data/atlas_branch/tlas_exome_bed_v1.4.3_allsamples_sm.vcf')
newvcf = open('ts_exomes_data/atlas_branch/atlas_exome_bed_v1.4.3_allsamples_sm_chrid.vcf', 'w')

head = oldvcf.next()
newvcf.write(head)
while head.startswith('##'):
    head = oldvcf.next()
    newvcf.write(head)
    

for line in oldvcf:
    line = 'chr' + line
    newvcf.write(line)

