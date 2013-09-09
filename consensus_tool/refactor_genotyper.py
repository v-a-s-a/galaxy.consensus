from ensemble_walker import *
from variant_ensemble import *

def __main__():

  testFiles = ['../data/atlas.small.vcf',
               '../data/freebayes.small.vcf',
               '../data/gatk.small.vcf']
  walker = concordant_walker(vcfList = testFiles)
  for x in walker.walk_concordant_sites():
    consensus = variant_ensemble(recordSet=x , samples=walker.samples)
    finalGenotypes = consensus.set_consensus()

if __name__ == '__main__':
  __main__()



