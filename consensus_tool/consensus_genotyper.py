#from ensemble_walker import concordant_walker
from vcf_ensemble import vcf_ensemble
from variant_ensemble import variant_ensemble
from consensus_writer import consensus_vcf
import argparse as arg
import pysam

def __main__():

  ## parse command line
  parser = arg.ArgumentParser(description='Find sites and genotypes which aggree among an arbitrary number of VCF files.')
  parser.add_argument('vcfFiles', nargs='+', metavar='VCFS', help='List of VCF files for input.')
  args = parser.parse_args()
  
  ## init the VCF ensemble
  ensemble = vcf_ensemble(vcfList = args.vcfFiles)
  
  outVcf = consensus_vcf()
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status")
  outVcf.add_format(id="GT", number="1", type="String", description="Genotype")
  outVcf = consensus_vcf()
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status")
  outVcf.add_format(id="GT", number="1", type="String", description="Genotype")
  outVcf.add_info(id="X1", number="1", type="String", description="Placeholder for INFO parsing")
  outVcf.add_info(id="X2", number="1", type="String", description="Placeholder 2 for INFO parsing")
  outVcf.samples = ensemble.samples
  outVcf.write_header()

  ## iterate over the concordant sites 
  for records, genotypes in ensemble.concordant_variants(siteThresh=3, genoThresh=3):
    outVcf.write_record( records, genotypes )

if __name__ == '__main__':
  __main__()



