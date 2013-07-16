
import sys
import optparse as opt
import sqlite3 as sql
import vcf as pyvcf
from pyVCF_record import _Record, _AltRecord, _Substitution
from consensus_functions import *
from genotyper import genotyper

def main():


  ## parse input arguments
  parser = opt.OptionParser()
  parser.add_option('--out', dest = 'baseOut', action = 'store', 
      help = 'File path base for output of .vcf and .log files.')
  parser.add_option('--atlas-vcf', dest = 'atlasVcf', action = 'store',
      help = 'Location of ATLAS vcf file for consensus.')
  parser.add_option('--gatk-vcf', dest = 'gatkVcf', action = 'store',
      help = 'Location of GATK vcf file for consensus.')
  parser.add_option('--freebayes-vcf', dest = 'freebayesVcf', action = 'store',
      help = 'Location of freebayes vcf file for consensus.')
  parser.add_option('--db-file', dest = 'dbFile', action = 'store',
      help = 'Location of file for sqlite db.')
  (options, args) = parser.parse_args()


  ## initialize the backend
  sqliteConnection = sql.connect(options.dbFile)

  ## make a consensus call set object
  consensusGenotyper = genotyper(sqliteConnection=sqliteConnection,
    atlasVcf=options.atlasVcf,
    gatkVcf=options.gatkVcf,
    freebayesVcf=options.freebayesVcf)

  ## create consensus genotype table
  consensusGenotyper.call_consensus(consThresh=3)

  ## write to file
  consensusGenotyper.write_vcf(options.baseOut)



if __name__ == '__main__':
  sys.exit(main())
