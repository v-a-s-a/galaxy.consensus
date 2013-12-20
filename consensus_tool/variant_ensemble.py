import itertools as it

class variant_ensemble:
  '''
  A set of variants from various VCF files. We want to be able to call various forms of consensus.

  A set will consist of pyVcf record objects.

  For now, we assume that all sites are concordant among all callers, i.e. we have a number of genotype sets equal to the number of call sets passed into consensus.

  '''   
  def __init__(self, *args, **kwargs):
    self.recordSet = kwargs.get('recordSet')
    self.samples = kwargs.get('samples') 
    self.ignoreMissing = kwargs.get('ignoreMissing')
    self.threshold = kwargs.get('threshold')

  def is_genotype_consensus(self, calls):
    '''
    Determine whether all _Calls for a sample at a site are equal.
    '''
    ## check to see if all calls are equal
    if calls[1:] == calls[:-1]: return True
    else: return False

  def concordant_genotypes(self):
    '''
    Return concordant genotypes at an arbitrary threshold.
    '''
    consensusGenotypes = dict()
    for sample in self.samples:
      genotypes = [ record.genotype(sample) for record in self.recordSet ]
      ## handle missing genotypes according to the flag passed in
      if self.ignoreMissing: genotypes = [ x for x in genotypes if x.gt_type != None ]
      if not genotypes:
        consensusGenotypes[sample] = './.'
        continue

      ## deliver the concordant genotypes


  def set_consensus(self):
    '''
    Return a consensus set of genotypes.
      Discount missing data during the vote.
    '''
    consensusGenotypes = dict()
    for sample in self.samples:
      genotypes = [ record.genotype(sample) for record in self.recordSet ]
      genotypes = [ x for x in genotypes if x.gt_type != None ]
      ## handle the missing genotype data here
      if not genotypes:
        consensusGenotypes[sample] = './.'
        continue
      consensusGenotype = self.is_genotype_consensus(calls=genotypes)
      if consensusGenotype: 
        ## return the consensus genotype
        consensusGenotypes[sample] = genotypes[0].gt_type
      else:
        ## no consensus so we return a '*' in its place
        consensusGenotypes[sample] = '*'
    return consensusGenotypes



