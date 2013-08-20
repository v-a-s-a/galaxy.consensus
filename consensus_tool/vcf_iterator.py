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
    
    print self.vcfs
   
    ## TODO: check that files are sorted
    for vcf in self.vcfs:
      pass


  def walk_concordant_sites(self):
    '''
    Return variants which match across call sets as you walk files
    '''

    inits = [ x.next() for x in self.readers ] 
    return inits


def __main__():
  
  testFiles = ['../data/atlas.small.vcf',
               '../data/freebayes.small.vcf',
               '../data/gatk.small.vcf']
  walker = concordant_walker(vcfList = testFiles)
  print walker.walk_concordant_sites()


if __name__ == '__main__':
  __main__()


