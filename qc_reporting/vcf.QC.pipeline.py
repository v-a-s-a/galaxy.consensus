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


## CLASSES
class vcfReport():
    '''
    A quality report generated from the vcf file.
    '''

    def __init__(self, vcfFile):
        self.vcfReader = vcf.Reader(vcfFile, strict_whitespace=True)
        #: number of transitions
        self.ts = 0
        #: list of MAF
        self.mafDist = list()
        #: number of missing
        #: list sample consensus-calling rate
        self.sampleConsenRate = list()
        #: list of site consensus-calling rates
        self.siteConsenRate = list()
        #: total number of called genotypes
        self.genotypesCalled = 0
        #: {sample ID : sample # of calls}
        self.sampleCallTotal = dict()
        for sample in self.vcfReader.samples:
            self.sampleCallTotal[sample] = 0
        
        ## go ahead an generate the report upon initialization
        self.__generateReport()


    def __generateReport(self):
        '''
        Iterate and process records in the vcf file.
        '''
        for record in self.vcfReader:
            self.__processRecord(record)

    def __processRecord(self, record):
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

        ## find samples with called genotypes
        #print self.sampleCallTotal
        #print record.num_called
        self.genotypesCalled += record.num_called
            
        ## is site a transition variant?
        if record.is_transition:
            self.ts += 1

        ## calculate and record maf
        if record.get_hom_refs() <= record.get_hom_alts():
            maf = record.aaf
        else:
            maf  = 1.0 - record.aaf
        self.mafDist.append(maf)

        ## calculate and record missingness per variant

        ## calculate and record mendelian errors for each call
        


## FUNCTIONS

## MAIN

def main():

    ## parse input arguments
    parser = opt.OptionParser()
    parser.add_option('--vcf', dest = 'vcfFile', help = 'Target VCF file.')
    (options, args) = parser.parse_args()
    vcfFile = options.vcfFile    

    ## create a QC report object
    qcreport = vcfReport(open(vcfFile, 'r'))

    ## TODO :: we want to wrap this in a nice PDF report
    ## potential libraries:
    ##  - sphinx
    ##  - pod
    ##  - reportLab

#    ## open connection to vcf
#    print vcfFile
#    reader = vcf.Reader(open(vcfFile), compressed = True)
#
#
#    for rec in reader:
#
#        ## process record. update aggregate metrics
#
#        print rec.get_unknowns()
#        print 'Call rate: %s' % rec.call_rate
#        print 'Unknown genotypes: %s' % rec.num_unknown
#        print 'Number of heterozygous calls: %i' % rec.num_het
#        print 'Total number called: %i' % rec.num_called
#        print 'Is transition: %r' % rec.is_transition

if __name__ == "__main__":
    main()
 






