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
  args = parser.parse_args()
 
  ## restrict current version to consensus-only thresholds 
  if args.genoThresh != len(args.vcfFiles):   
    print "STOPPING"
    print "\tGenotype threshold does not match the number of input VCF files. Non-consensus modes of ensemble calling are not yet supported."
    sys.exit()
  if args.siteThresh != len(args.vcfFiles):
    print "STOPPING"
    print "\tVariant site threshold does not match the number of input VCF files. Non-consensus modes of ensemble calling are not yet supported."
    sys.exit()

 
  ## create the VCF ensemble
  ensemble = vcf_ensemble(vcfList = args.vcfFiles)
 
  ## setup output VCF file. Dummy fields are created for downstream parsing with other tools.
  outVcf = consensus_vcf()
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status")
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



