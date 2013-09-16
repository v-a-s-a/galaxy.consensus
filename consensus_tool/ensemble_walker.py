'''
Given a list of vcf files, iterate over records in each simultaneously.

** refactor this into a class with an __iter__ method w/ yields
** use the tabix indexing to access regions of files

'''

import vcf as pyvcf

class concordant_walker:
  '''
  Walks across an arbitrary number of VCF files given a chomosome or contig contained in all
  '''


  @property
  def samples(self):
    return self.samples

  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    self.contig = kwargs.get('contig')
    self.readers = [ pyvcf.Reader(open(x), prepend_chr=True) for x in self.vcfs ]
    self.readers = [ x.fetch(chrom=self.contig, start=1, end=1000000000) for x in self.readers ]

    ## number of readers passed to this instance
    self.vcfCount = len(self.vcfs)
    ## first records in each file
    self.records = [ x.next() for x in self.readers ] 
    ## track current positions of each record
    self.positions = [ int(rec.POS) for rec in self.records  ]
    ## track current greatest position (known as the prime position, and prime readers)
    self.prime = max(self.positions)
    ## pull inidices of readers which have reached the largest observed position
    self.primeIndices = [ i for i in range(self.vcfCount) if self.positions[i]==self.prime ]
    ## identify samples which match between sets
    sampleSets = [ set(x.samples) for x in self.readers  ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )

    ## TODO: check that files are sorted
    for vcf in self.vcfs:
      pass

  def __updatePrimeIndices(self):
    '''
    Update indices of VCF readers which match the current largest position.
    '''
    self.primeIndices = [ i for i, x in enumerate(self.positions) if x==self.prime ]


  def __updateRecords(self):
    '''
    Iterate non-prime readers.
    '''
    ## pull indices of readers which need to iterate to reach or exceed the prime position 
    iterIndices = [ i for i in range(self.vcfCount) if i not in self.primeIndices ] 
    ## update records in non-prime readers
    self.records = [ x if i in self.primeIndices else self.readers[i].next() for i,x in enumerate(self.records) ]
    self.positions = [ int(x.POS) for x in self.records ]
    ## update prime if necessary
    if max(self.positions) > self.prime:
      self.prime = max(self.positions)
    ## update which readers we are iterating over
    self.__updatePrimeIndices()


  def __getNewPositions(self):
    '''
    Try a new set of positions after processing concordant sites.
    '''
    self.positions = [ int(x.next().POS) for x in self.readers ]
    self.__updatePrimeIndices()

  def __iterUntilEqual(self):
    '''
    Update the positions vector until all positions are equal.
    '''
    ## check if all readers have reached the same position
    while len(set(self.positions)) != 1:
      try:
        self.__updateRecords()
      except StopIteration:
        ## we've exhausted one of the readers
        return False
    ## we have not yet exhausted the readers
    return True

  def walk_concordant_sites(self):
    '''
    Return variants which match across call sets as you walk files
    
    Start simple:
      * assume vcf is sorted
      * assume only SNPs
      * assume only bi-allelic variation
      * assume only one chromosome
    '''

    while self.__iterUntilEqual():
      ## move on to next set of variants
      self.__getNewPositions()

      ## yield set genotypes
      yield self.records

    


