import vcf as pyvcf
import itertools as it
from variant_ensemble import variant_ensemble

class simple_variant:
  '''
  Minimum representation of a variant.
  '''

  def __init__(self, rec, ALT):
    '''
    Take a pyVCF _Record and make a string. Must provide explicit alternate allelel ALT.
    '''
    self.ID = rec.CHROM + ':' + str(rec.POS) + ':' + str(rec.REF) + ':' + str(ALT)
    self.contig = rec.CHROM
    self.pos = rec.POS
  def __hash__(self):
    return hash(self.ID)
  def __eq__(self, other):
    if self.ID == other.ID: return True
    else: return False
  def __str__(self):
    return self.ID    


class vcf_ensemble:
  '''
  Represents an arbitrary number of VCF files describing the same data set.
  '''
  @property
  def samples(self):
    '''
    Samples common to all input VCF files.
    '''
    sampleSets = [ set(x.samples) for x in self.vcfReaders ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )
    return self.samples

  @samples.setter
  def samples(self, samples):
    self.samples = samples

  @samples.getter
  def samples(self):
    sampleSets = [ set(x.samples) for x in self.vcfReaders ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )
    return self.samples
  
  @property
  def variants(self):
    '''
    List of variant sets found in each caller.
    '''
    variantSets = []
    readers = [ pyvcf.Reader(open(x), prepend_chr = False) for x in self.vcfs]
    for x in readers:
      variants = set()
      for rec in x:
        for alt in rec.ALT:
          variants.add(simple_variant(rec, alt))

      variantSets.append(variants)
    return variantSets

  @property
  def vcfReaders(self):
    '''
    Connections to VCF files.
    '''
    return [ pyvcf.Reader(open(x), prepend_chr = False) for x in self.vcfs ]


  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    self.ignoreMissing = kwargs.get('ignoreMissing')

  def _concordant_sites(self, threshold):
    '''
    Find variants that agree at some threshold.

    This problem reduces to finding the union of a set of sets. The top level
    set consists of combinations of VCF sets set by the threshold (i.e.
   variants choose threshold).
    '''
    sets = set()
    for comb in it.combinations(self.variants, threshold):
      sets.add( frozenset(reduce( lambda x, y: x.intersection(y), comb )) ) 
    return reduce( lambda x, y: x.union(y), sets )

  def concordant_variants(self, siteThresh, genoThresh):
    '''
    Reduce the ensemble of VCFs into a set of variant sites that agree at some threshold.
    '''
    variants = self._concordant_sites(siteThresh)
    samples = self.samples
    for variant in variants:
      records = [ x.fetch(variant.contig, variant.pos) for x in self.vcfReaders ]
      ## clean up missing records for a given threshold
      records = [ x for x in records if x ]
      ensemble = variant_ensemble(recordSet=[ x for x in records if x], samples=samples, threshold = genoThresh, ignoreMissing = self.ignoreMissing)
      yield records, ensemble.set_concordance()



