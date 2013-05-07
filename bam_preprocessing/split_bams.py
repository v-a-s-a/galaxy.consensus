'''
Split an entire sequencing experiment into chromosomes.

Input:
  Directory of sample BAM files

Output:
  A directory of directories 23 + X + Y containing subsetted BAM files for each

Dependencies:
  bamtools?
  samtools?

TODO:
  * parse inputs from command line.
  * thread or otherwise parallelize this.

'''

import os
import subprocess as sp

def __main__():

  ## process command line inputs
  bamDir = '/path/to/experiment/bams/'
  metaDir = '/path/to/chromosome/directories/'
  bamtoolsLoc = '/path/to/bamtools/binary'


  ## store BAM file paths
  for root, dirs, files in os.walk(bamDir):
    bamPaths = [root + f for f in files]

  ## create chromosome directories
  chroms = range(1,23).append('X').append('Y')
  chroms = [ str(field) for field in chroms ]
  for chrom in chroms:
    chromDir = metaDir + 'chrm_' + chrom + '/'
    if not os.path.exists(chromDir):
      os.makedirs(chromDir)
  
  ## split each BAM file
  for bam in bamPaths:
    splitCall = [bamtoolsLoc, '-in', '-reference']
    sp.Call(splitCall)
    ## move chromosome BAMs into appropriate directory
    for chrom in chroms:
      mvCall = ['mv']


if __name__ == '__main__':
  __main__()

