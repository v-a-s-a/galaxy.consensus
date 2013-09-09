class consensus_writer:
  '''
  Take [site info, dicts of sample:genotype] and write out a valid VCF file.

  Currently, the pyVcf module wants you to provide a template reader to write. This is not a very good solution for a situation in which you are writing a brand new VCF file.

  Additionally, I would like to easily pass data from input VCF files and have them written to a specified field in the final VCF file.
  '''

class consensus_vcf:

  @property
  header

  @property
  format

  @property
  info


