import vcf as pyvcf
import pysam

def get_regions(vcf):
  '''
  Function to find (contig, start, end) for contig in VCF
  '''

  reader = pyvcf.Reader(open(vcf))
  tabixReader = pysam.Tabix(vcf)

  regions = list()
  for contig in tabixReader.contigs:
    row = tabixReader._parseRegion(region=contig)
    start = row[-2]
    end = row[-1]
    print (contig, start, end)
  
  



