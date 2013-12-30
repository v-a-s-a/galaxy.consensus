'''
Given a list of vcf files, iterate over records in each simultaneously.

** refactor this into a class with an __iter__ method w/ yields
** use the tabix indexing to access regions of files

'''

import vcf as pyvcf

class concordant_walker:
  '''
  Walks across an arbitrary number of VCF files given a common contig.
  '''

  @property
  def positions(self):
    '''
    Current positions in each VCF file.
    '''
    return [ x.POS for x in self.records ]

  @property
  def primeIndices(self):
    '''
    Indices of VCF readers which match the current largest position.
    '''
    return [ i for i, x in enumerate(self.positions) if x==self.prime ]

  @property
  def prime(self):
    '''
    Largest position currently observed.
    '''
    return max(self.positions)

  @property
  def vcfCount(self):
    '''
    Number of VCF files passed in to this instance of the consensus walker.
    '''
    return len(self.vcfs)

  @property
  def samples(self):
    '''
    Samples common to all input VCF files.
    '''
    sampleSets = [ set(x.samples) for x in self.readers  ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )
    return self.samples
  

  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    ## start readers across the VCF files
    self.readers = [ pyvcf.Reader(open(x), prepend_chr = False) for x in self.vcfs ]
    
    self._getNewRecords()
    ## TODO: check that files are sorted
    for vcf in self.vcfs:
      pass

  def _updateRecords(self):
    '''
    Iterate non-prime readers.
    '''
    ## update records in non-prime readers
    self.records = [ x if i in self.primeIndices else self.readers[i].next() for i,x in enumerate(self.records) ]
    ## update prime if necessary
    if max(self.positions) > self.prime:
      self.prime = max(self.positions)

  def _getNewRecords(self):
    '''
    Get a new set of records after processing concordant sites or initilization.
    '''
    self.records = [ x.next() for x in self.readers ]
    self._updateRecords()


  def walk_concordant_sites(self):
    '''
    Return variants which match across call sets as you walk files.
    
    Start simple:
      * assume only SNPs
      * assume only bi-allelic variation
    '''
    ## iterate until all positions are equal
    while True:
      if len(set(self.positions)) != 1:
        try: self._updateRecords()
        ## stop only when a reader has been exhausted
        except StopIteration:
          break
      else:
        yield self.records
        self._getNewRecords()
 

