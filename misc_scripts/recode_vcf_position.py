'''
recode the chromosomal identifier in a vcf
'''

oldvcf = open('ts_exomes_data/atlas_exome.vcf')
newvcf = open('ts_exomes_data/atlas_exome_chrid.vcf', 'w')

head = oldvcf.next()
newvcf.write(head)
while head.startswith('##'):
    head = oldvcf.next()
    newvcf.write(head)
    

for line in oldvcf:
    line = 'chr' + line
    newvcf.write(line)

