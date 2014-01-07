#from ensemble_walker import concordant_walker
from vcf_ensemble import vcf_ensemble
from variant_ensemble import variant_ensemble
from consensus_writer import consensus_vcf
import argparse as arg
import pysam
import sys

def __main__():

  ## parse command line
  parser = arg.ArgumentParser(description='Find sites and genotypes that aggree among an arbitrary number of VCF files.')
  parser.add_argument('vcfFiles', nargs='+', metavar='VCFS', help='List of VCF files for input.')
  parser.add_argument('--site-threshold', '-s', dest='siteThresh', action='store', type=int, help='Number of inputs which must agree for a site to be included in the output.')
  parser.add_argument('--genotype-threshold', '-g', dest='genoThresh', action='store', type=int, help='Number of inputs which must agree for a genotype to be marked as non-missing.')
  parser.add_argument('--ignore-missing', '-m', dest='ignoreMissing', action='store_true', help='Flag specifying how to handle missing genotypes in the vote. If present, missing genotypes are excluded from the genotype concordance vote unless all genotypes are missing.')
  args = parser.parse_args()
 
  ## create the VCF ensemble
  ensemble = vcf_ensemble(vcfList = args.vcfFiles, ignoreMissing = args.ignoreMissing)
 
  ## setup output VCF file. Dummy fields are created for downstream parsing with other tools.
  outVcf = consensus_vcf()
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status. \'C\' is concordant, \'D\' is discordant, and \'A\' is ambiguous (no clear agreement).")
  outVcf.add_format(id="GT", number="1", type="String", description="Genotype")
  outVcf.add_info(id="X1", number="1", type="String", description="Placeholder for INFO parsing")
  outVcf.add_info(id="X2", number="1", type="String", description="Placeholder 2 for INFO parsing")
  outVcf.samples = ensemble.samples
  outVcf.write_header()

  ## iterate over the concordant sites 
  for records, genotypes in ensemble.concordant_variants(siteThresh=args.siteThresh, genoThresh=args.genoThresh):
    outVcf.write_record( records, genotypes )


if __name__ == '__main__':
  __main__()



