'''
Given a list of vcf files, iterate over records in each simultaneously.

Check that files are sorted.




** refactor this into a class with an __iter__ method w/ yields
** use the tabix indexing to access regions of files

'''

import vcf as pyvcf

class concordant_walker:

  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    self.readers = [ pyvcf.Reader(open(x)) for x in self.vcfs ]
    ## first records in each file
    self.inits = [ x.next() for x in self.readers ] 
    ## track current positions of each record
    self.positions = [ int(rec.POS) for rec in self.inits  ]
    ## track current greatest position
    self.prime = max(self.positions)
   
    ## TODO: check that files are sorted
    for vcf in self.vcfs:
      pass

    def __updatePositions(self):
      '''
      Update the positions vector until all positions are equal, or
      a new prime position is found.
      '''
      
      
      for reader in self.readers:
        


  def walk_concordant_sites(self):
    '''
    Return variants which match across call sets as you walk files
    
    Start simple:
      * assume vcf is sorted
      * assume only SNPs
      * assume only bi-allelic variation
      * assume only one chromosome

    '''
    ## find largest site among these, X
    ## iterate other files until records match or are exceeded
    ## if records match -- do consensus work/exit/whatever
    ## if record found exceeding position of X, set this as X and iterate

    return positions, prime


def __main__():
  
  testFiles = ['../data/atlas.small.vcf',
               '../data/freebayes.small.vcf',
               '../data/gatk.small.vcf']
  walker = concordant_walker(vcfList = testFiles)
  print walker.walk_concordant_sites()


if __name__ == '__main__':
  __main__()


