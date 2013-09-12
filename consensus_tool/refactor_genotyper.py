from ensemble_walker import concordant_walker
from variant_ensemble import variant_ensemble
from consensus_writer import consensus_vcf
import argparse as arg

def __main__():

  ## parse command line
  parser = arg.ArgumentParser(description='Find sites and genotypes which aggree among an arbitrary number of VCF files.')
  parser.add_argument('vcfFiles', nargs='+', metavar='VCFS', help='List of VCF files for input.')
  args = parser.parse_args()

  ## instantiate a walker over the input vcf files.
  walker = concordant_walker(vcfList = args.vcfFiles)

  ## set up the consensus VCF you want to write out
  ## TODO:: there should be a standard and transparent way to propagate information for individual VCF files to the consensus stage.
  outVcf = consensus_vcf()
  ## leave these for later -- there is a better way
  #outVcf.add_info(id="AQ", number="1", type="Float", description="Atlas QUAL score")
  #outVcf.add_info(id="FQ", number="1", type="Float", description="Freebayes QUAL score")
  #outVcf.add_info(id="GQ", number="1", type="Float", description="GATK QUAL score")
  #outVcf.add_format(id="PL", number="G", type="Integer", description="GATK's normalized, Phred-scaled likelihoods for genotypes as defined in their VCF spec")
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status")
  outVcf.add_format(id="GT", number="1", type="String", description="Genotype")
  outVcf.samples = walker.samples
  outVcf.write_header()


  ## iterate over sites which match among all VCFs in concordant_walker
  for x in walker.walk_concordant_sites():
    consensus = variant_ensemble(recordSet=x , samples=walker.samples)
    finalGenotypes = consensus.set_consensus()
    outVcf.write_record(x, finalGenotypes)


if __name__ == '__main__':
  __main__()



