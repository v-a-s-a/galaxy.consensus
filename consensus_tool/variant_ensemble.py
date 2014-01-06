import itertools as it
from collections import Counter
from my_exceptions import ambiguousConcordance
from my_exceptions import discordant

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

  def _genotype_concordance(self, calls, thresold):
    '''
    Find the call which agrees at a certain threshold. If a tie is observed, an exception is raised, and a missing value will be written to the VCF file and flagged as ambiguous.
    '''
    ## count occurences of a genotype, and those which match the threshold
    counted = Counter(calls).most_common()
    matched = [ x[0] for x in counted if x[1]>=self.threshold ]
    ## check for a tie -- and raise an exception if necessary
    if len(matched) > 1: raise ambiguousConcordance('tie!')
    elif not matched: raise discordant('no match')
    else: return matched[0]
    
  def set_concordance(self, threshold):
    '''
    Return a set of concordant genotypes at a specified threshold.
    '''

    concordantGenotypes = dict()
    for sample in self.samples:
      calls = [ record.genotype(sample).gt_type for record in self.recordSet ]
      #genotypes = [ x for x in genotypes if x.gt_type != None ] ## now treating missing as valid genotype
      ## handle the missing genotype data here
      if not calls:
        concordantGenotypes[sample] = './.'
        continue

      ## try to find a concordant genotype
      try:
        concordantGenotypes[sample] = self._genotype_concordance(calls, threshold)
      except ambiguousConcordance:
        ## really we should allow this to bubble up and be handled when writing the VCF
        concordantGenotypes[sample] = '**'
      except discordant:
        concordantGenotypes[sample] = '*'
     

    return concordantGenotypes

