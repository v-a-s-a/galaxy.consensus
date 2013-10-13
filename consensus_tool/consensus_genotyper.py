from ensemble_walker import concordant_walker
from variant_ensemble import variant_ensemble
from consensus_writer import consensus_vcf
import argparse as arg

def __main__():

  ## parse command line
  parser = arg.ArgumentParser(description='Find sites and genotypes which aggree among an arbitrary number of VCF files.')
  parser.add_argument('vcfFiles', nargs='+', metavar='VCFS', help='List of VCF files for input.')
  parser.add_argument('--contig', help='Contig ID to operate on.', required = False)
  args = parser.parse_args()

  print args.contig

<<<<<<< HEAD
  if args.contig:
    contigs = [args.contig]  
  else:
    contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] 

  walker = concordant_walker(vcfList = args.vcfFiles, contig = contigs[0])
=======
  walker = concordant_walker(vcfList = args.vcfFiles, contig = '22')
>>>>>>> e38edd9045d8ed7519d05641fb8aa5e894086bd6
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
  outVcf.add_info(id="X1", number="1", type="String", description="Placeholder for INFO parsing")
  outVcf.add_info(id="X2", number="1", type="String", description="Placeholder 2 for INFO parsing")
  outVcf.samples = walker.samples
  outVcf.write_header()


  ## walk over each contig independently
<<<<<<< HEAD
  #contigs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrMT']
=======
  #contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] 
  contigs = ['chr22']
>>>>>>> e38edd9045d8ed7519d05641fb8aa5e894086bd6

  for contig in contigs:
    ## instantiate a walker over the input vcf files.
    walker = concordant_walker(vcfList = args.vcfFiles, contig = contig)

    ## iterate over sites which match among all VCFs in concordant_walker
    for concordantSites in walker.walk_concordant_sites():
      consensus = variant_ensemble(recordSet=concordantSites, samples=walker.samples)
      concordantGenotypes = consensus.set_consensus()
      outVcf.write_record(concordantSites, concordantGenotypes)


if __name__ == '__main__':
  __main__()



