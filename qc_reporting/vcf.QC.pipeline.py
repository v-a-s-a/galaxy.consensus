'''
A pipeline for generating a set of QC metrics for an arbitrary VCF.

Target Metrics:
    Current:
        * ts/tv ratio
        * sample missingness rate
        * site missingness rate
        * sample non-consensus rate
        * site non-consensus rate
        * MAF distribution

    Planned:
        * combined QUAL
        * mendelian error

programs used:
    * vcftools
    * PLINK


'''

## MODULES
import vcf
import optparse as opt
import subprocess as sp
import numpy as np
import sys

## CLASSES
class vcfReport():
    '''
    A quality report generated from the vcf file.
    '''

    def __init__(self, vcfFile, mvar):
        self.vcfReader = vcf.Reader(open(vcfFile, 'r'))
        self.nsam = len(self.vcfReader.samples)
        self.mvar = mvar
        self.vcfFile = vcfFile
        #: number of transitions
        self.ts = 0
        #: number of transversions
        self.tv = 0
        #: list of MAF
        self.mafDist = np.empty( self.mvar, dtype=float)
        #: list of call rates
        self.callrateDist = np.empty( self.mvar, dtype=float)
        #: list of site consensus-calling rates
        self.siteConsRate = np.empty( self.mvar, dtype=float)
        #: total number of called genotypes
        self.genotypesCalled = 0
        #: {sample ID : sample # of calls}
        self.sampleCallTotal = np.zeros( (self.nsam, 1), dtype=int)
        
        ## go ahead an generate the report upon initialization
        self.__generateReport()


    def __generateReport(self):
        '''
        Iterate and process records in the vcf file.
        Calculate mendelian errors using PLINK.
        '''
        for i, record in enumerate(self.vcfReader):
            self.__processRecord(record, i)

        self.__findMendelErrors()

    def findMendelErrors(self):
      '''
      Reformat VCF as plink files, and add pedigree information.
      Calculate the rates of mendelian errors from there.
      
      Inputs:
          - output directory for intermediate files
          - file mapping SM to RG for our samples
          - file containing pedigree information

      '''
      


    def __processRecord(self, record, recidx):
        '''
        Update sample level statistics:
            * total number of called genotypes
            * non-consensus rate
            * mendelian error rate

        Store site level statistics:
            * missingness rate
            * non-consensus rate
            * mendelian error rate
            * minor allele frequency
        '''


        ## per sample totals
        consCount = 0
        for idx, sample in enumerate(record.samples):
          ## check if genotype is called
          if sample.gt_type:
              self.sampleCallTotal[idx] += 1
          ## check if genotype was called consensus
          if getattr(sample.data ,'CN') == 'T':
              consCount += 1

        ## global call total
        self.genotypesCalled += record.num_called
            
        ## is site a transition variant?
        if record.is_transition:
            self.ts += 1
        else:
            self.tv += 1

        ## calculate and record maf
        if record.get_hom_refs() <= record.get_hom_alts():
            maf = record.aaf
        else:
            maf  = 1.0 - record.aaf
        self.mafDist[recidx] = maf

        ## calculate and record missingness per variant
        self.callrateDist[recidx] = record.call_rate

        ## if consensus rate is available -- record a consensus rate
        self.siteConsRate[recidx] = consCount / record.num_called
        

### FUNCTIONS ###
#################

## quickly find the number of samples in a vcf
def num_variants():
    ## grep -v '^#' test.consensus.vcf | wc -l
    return Null 

## quickly find the number of variants in a vcf
def num_samples(vcfReader):
    return len(vcfReader.samples)


### MAIN ###
############

def main():

    ## parse input arguments
    parser = opt.OptionParser()
    parser.add_option('--vcf', dest = 'vcfFile', help = 'Target VCF file.')
    parser.add_option('--num-variants', dest = 'mvar', help = 'Number of variants in target VCF file.')
    (options, args) = parser.parse_args()
    vcfFile = options.vcfFile    
    mvar = int(options.mvar)

    ## create a QC report object
    qcreport = vcfReport(vcfFile, mvar)

    ## TODO :: we want to wrap this in a nice PDF report
    ## potential libraries:
    ##  - sphinx
    ##  - pod
    ##  - reportLab


if __name__ == "__main__":
    main()
 






