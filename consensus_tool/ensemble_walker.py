'''
Given a list of vcf files, iterate over records in each simultaneously.

** refactor this into a class with an __iter__ method w/ yields
** use the tabix indexing to access regions of files

'''

import vcf as pyvcf

class concordant_walker:

  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    self.readers = [ pyvcf.Reader(open(x)) for x in self.vcfs ]
    ## number of readers passed to this instance
    self.vcfCount = len(self.vcfs)
    ## first records in each file
    self.inits = [ x.next() for x in self.readers ] 
    ## track current positions of each record
    self.positions = [ int(rec.POS) for rec in self.inits  ]
    ## track current greatest position (known as the prime position, and prime readers)
    self.prime = max(self.positions)
    ## pull inidices of readers which have reached the largest observed position
    self.primeIndices = [ i for i in range(self.vcfCount) if self.positions[i]==self.prime ]

    ## TODO: check that files are sorted
    for vcf in self.vcfs:
      pass

  def __updatePrimeIndices(self):
    '''
    Update indices of VCF readers which match the current largest position.
    '''
    self.primeIndices = [ i for i, x in enumerate(self.positions) if x==self.prime ]


  def __updatePositions(self):
    '''
    Iterate non-prime readers.
    '''
    ## pull indices of readers which need to iterate to reach or exceed the prime position 
    iterIndices = [ i for i in range(self.vcfCount) if i not in self.primeIndices ] 
    ## update positions in non-prime readers
    self.positions = [ x if i in self.primeIndices else int(self.readers[i].next().POS) for i,x in enumerate(self.positions) ]
    ## update prime if necessary
    if max(self.positions) > self.prime:
      self.prime = max(self.positions)
    ## update which readers we are iterating over
    self.__updatePrimeIndices()



    print self.positions, self.primeIndices

  def __iterUntilEqual(self):
    '''
    Update the positions vector until all positions are equal.
    '''
    ## check if all readers have reached the same position
    while len(set(self.positions)) != 1:
      self.__updatePositions()


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

    self.__iterUntilEqual()
    return self.positions


def __main__():
  
  testFiles = ['../data/atlas.small.vcf',
               '../data/freebayes.small.vcf',
               '../data/gatk.small.vcf']
  walker = concordant_walker(vcfList = testFiles)
  print walker.walk_concordant_sites()


if __name__ == '__main__':
  __main__()


