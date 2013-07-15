'''
A class for a branch feeding into the consensus caller.

  * keep track of the source VCF file
  * interface with the sqlite backend
'''


class branch:
  
  def branch_init(vcfFile, tableName, dbCon):
    '''
    Refactor the create_vcf_tale function using a vcf parser

    Creates a table in a given sqlite db.

    Parameters:
        @vcfFile: a valid vcf for input into table
        @dbCon: connection to consensus sqlite db
        @tableName: name corresponding to vcf file

    Return: None

   db table format:
        -row is variant
        -col for varID s chr_pos
        -col for chromosome
        -col for position
        -col for REF allele
        -col for ALT allele
        -col for QUAL
        -col for sample genotype -- sample IDs are enclosed by "
    '''

  
