'''
Consensus caller for vcf files

the vcf structure is such that the first field in the file name delimited by
'.' is the name of the table in the db

'''

###############
### MODULES ###
###############

import sys
import sqlite3 as sql
import os
import optparse as opt
import vcf as vcfMod


###############
### CLASSES ###
###############


## here I stole an internal class from james casbon's pcvcf
class _Record(object):
    """ A set of calls at a site.  Equivalent to a row in a VCF file.

        The standard VCF fields CHROM, POS, ID, REF, ALT, QUAL, FILTER,
        INFO and FORMAT are available as properties.

        The list of genotype calls is in the ``samples`` property.
    """
    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
            sample_indexes, samples=None):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        #: 0-based start coordinate
        self.start = self.POS - 1
        #: 1-based end coordinate
        self.end = self.start + len(self.REF)
        #: list of alleles. [0] = REF, [1:] = ALTS
        self.alleles = [self.REF]
        self.alleles.extend(self.ALT)
        #: list of ``_Calls`` for each sample ordered as in source VCF
        self.samples = samples
        self._sample_indexes = sample_indexes

    def __eq__(self, other):
        """ _Records are equal if they describe the same variant (same
        position, alleles) """
        return (self.CHROM == other.CHROM and
                self.POS == other.POS and
                self.REF == other.REF and
                self.ALT == other.ALT)

def __iter__(self):
        return iter(self.samples)

    def __str__(self):
        return "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s)" % self.__dict__

    def __cmp__(self, other):
        return cmp((self.CHROM, self.POS), (other.CHROM, other.POS))

    def add_format(self, fmt):
        self.FORMAT = self.FORMAT + ':' + fmt

    def add_filter(self, flt):
        self.FILTER.append(flt)

    def add_info(self, info, value=True):
        self.INFO[info] = value

    def genotype(self, name):
        """ Lookup a ``_Call`` for the sample given in ``name`` """
        return self.samples[self._sample_indexes[name]]

    @property
    def num_called(self):
        """ The number of called samples"""
        return sum(s.called for s in self.samples)

    @property
    def call_rate(self):
        """ The fraction of genotypes that were actually called. """
        return float(self.num_called) / float(len(self.samples))

    @property
    def num_hom_ref(self):
        """ The number of homozygous for ref allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 0])

    @property
    def num_hom_alt(self):
        """ The number of homozygous for alt allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 2])

    @property
    def num_het(self):
        """ The number of heterozygous genotypes"""
        return len([s for s in self.samples if s.gt_type == 1])

    @property
    def num_unknown(self):
        """ The number of unknown genotypes"""
        return len([s for s in self.samples if s.gt_type is None])

    @property
    def aaf(self):
        """ The allele frequency of the alternate allele.
           NOTE 1: Punt if more than one alternate allele.
           NOTE 2: Denominator calc'ed from _called_ genotypes.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        het = self.num_het
        hom_alt = self.num_hom_alt
        num_chroms = float(2.0 * self.num_called)
        return float(het + 2 * hom_alt) / float(num_chroms)

    @property
    def nucl_diversity(self):
        """
        pi_hat (estimation of nucleotide diversity) for the site.
        This metric can be summed across multiple sites to compute regional
        nucleotide diversity estimates.  For example, pi_hat for all variants
        in a given gene.

        Derived from:
        \"Population Genetics: A Concise Guide, 2nd ed., p.45\"
          John Gillespie.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        p = self.aaf
        q = 1.0 - p
        num_chroms = float(2.0 * self.num_called)
        return float(num_chroms / (num_chroms - 1.0)) * (2.0 * p * q)

    def get_hom_refs(self):
        """ The list of hom ref genotypes"""
        return [s for s in self.samples if s.gt_type == 0]

    def get_hom_alts(self):
        """ The list of hom alt genotypes"""
        return [s for s in self.samples if s.gt_type == 2]

    def get_hets(self):
        """ The list of het genotypes"""
        return [s for s in self.samples if s.gt_type == 1]

    def get_unknowns(self):
        """ The list of unknown genotypes"""
        return [s for s in self.samples if s.gt_type is None]

    @property
    def is_snp(self):
        """ Return whether or not the variant is a SNP """
        if len(self.REF) > 1:
            return False
        for alt in self.ALT:
            if alt is None or alt.type != "SNV":
                return False
            if alt not in ['A', 'C', 'G', 'T']:
                return False
        return True

    @property
    def is_indel(self):
        """ Return whether or not the variant is an INDEL """
        is_sv = self.is_sv

        if len(self.REF) > 1 and not is_sv:
            return True
        for alt in self.ALT:
            if alt is None:
                return True
            if alt.type != "SNV" and alt.type != "MNV":
                return False
            elif len(alt) != len(self.REF):
                # the diff. b/w INDELs and SVs can be murky.
                if not is_sv:
                    # 1 2827693 .   CCCCTCGCA   C   .   PASS    AC=10;
                    return True
                else:
                    # 1 2827693 .   CCCCTCGCA   C   .   PASS    SVTYPE=DEL;
                    return False
        return False

    @property
    def is_sv(self):
        """ Return whether or not the variant is a structural variant """
        if self.INFO.get('SVTYPE') is None:
            return False
        return True

    @property
    def is_transition(self):
        """ Return whether or not the SNP is a transition """
        # if multiple alts, it is unclear if we have a transition
        if len(self.ALT) > 1:
            return False

        if self.is_snp:
            # just one alt allele
            alt_allele = self.ALT[0]
            if ((self.REF == "A" and alt_allele == "G") or
                (self.REF == "G" and alt_allele == "A") or
                (self.REF == "C" and alt_allele == "T") or
                (self.REF == "T" and alt_allele == "C")):
                return True
            else:
                return False
        else:
            return False

    @property
    def is_deletion(self):
        """ Return whether or not the INDEL is a deletion """
        # if multiple alts, it is unclear if we have a transition
        if len(self.ALT) > 1:
            return False

        if self.is_indel:
            # just one alt allele
            alt_allele = self.ALT[0]
            if alt_allele is None:
                return True
            if len(self.REF) > len(alt_allele):
                return True
            else:
                return False
        else:
            return False

    @property
    def var_type(self):
        """
        Return the type of variant [snp, indel, unknown]
        TO DO: support SVs
        """
        if self.is_snp:
            return "snp"
        elif self.is_indel:
            return "indel"
        elif self.is_sv:
            return "sv"
        else:
            return "unknown"

@property
    def var_subtype(self):
        """
        Return the subtype of variant.
        - For SNPs and INDELs, yeild one of: [ts, tv, ins, del]
        - For SVs yield either "complex" or the SV type defined
          in the ALT fields (removing the brackets).
          E.g.:
               <DEL>       -> DEL
               <INS:ME:L1> -> INS:ME:L1
               <DUP>       -> DUP

        The logic is meant to follow the rules outlined in the following
        paragraph at:

        http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
        "For precisely known variants, the REF and ALT fields should contain
        the full sequences for the alleles, following the usual VCF conventions.
        For imprecise variants, the REF field may contain a single base and the
        ALT fields should contain symbolic alleles (e.g. <ID>), described in more
        detail below. Imprecise variants should also be marked by the presence
        of an IMPRECISE flag in the INFO field."
        """
        if self.is_snp:
            if self.is_transition:
                return "ts"
            elif len(self.ALT) == 1:
                return "tv"
            else:  # multiple ALT alleles.  unclear
                return "unknown"
        elif self.is_indel:
            if self.is_deletion:
                return "del"
            elif len(self.ALT) == 1:
                return "ins"
            else:  # multiple ALT alleles.  unclear
                return "unknown"
        elif self.is_sv:
            if self.INFO['SVTYPE'] == "BND":
                return "complex"
            elif self.is_sv_precise:
                return self.INFO['SVTYPE']
            else:
                return self.ALT[0].type
        else:
            return "unknown"

    @property
    def sv_end(self):
        """ Return the end position for the SV """
        if self.is_sv:
            return self.INFO['END']
        return None


@property
    def is_sv_precise(self):
        """ Return whether the SV cordinates are mapped
            to 1 b.p. resolution.
        """
        if self.INFO.get('IMPRECISE') is None and not self.is_sv:
            return False
        elif self.INFO.get('IMPRECISE') is not None and self.is_sv:
            return False
        elif self.INFO.get('IMPRECISE') is None and self.is_sv:
            return True

    @property
    def is_monomorphic(self):
        """ Return True for reference calls """
        return len(self.ALT) == 1 and self.ALT[0] is None


class _AltRecord(object):
    '''An alternative allele record: either replacement string, SV placeholder, or breakend'''
    __metaclass__ = ABCMeta

    def __init__(self, type, **kwargs):
        super(_AltRecord, self).__init__(**kwargs)
        #: String to describe the type of variant, by default "SNV" or "MNV", but can be extended to any of the types described in the ALT lines of the header (e.g. "DUP", "DEL", "INS"...)
        self.type = type

    @abstractmethod
    def __str__(self):
        raise NotImplementedError

    def __eq__(self, other):
        return self.type == other.type


class _Substitution(_AltRecord):
    '''A basic ALT record, where a REF sequence is replaced by an ALT sequence'''

    def __init__(self, nucleotides, **kwargs):
        if len(nucleotides) == 1:
            super(_Substitution, self).__init__(type="SNV", **kwargs)
        else:
            super(_Substitution, self).__init__(type="MNV", **kwargs)
        #: Alternate sequence
        self.sequence = str(nucleotides)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.sequence)

    def __eq__(self, other):
        if isinstance(other, basestring):
            return self.sequence == other
        else:
            return super(_Substitution, self).__eq__(other) and self.sequence
== other.sequence




#################
### FUNCTIONS ###
#################

def store_vcf(vcfFile, tableName, dbCon):
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
        -col for sample genotype -- sample IDs are enclosed by "
    '''

    ## initialize vcf parser
    reader = vcfMod.Reader(open(vcfFile, 'r'))


    ## form columns for addition into table
    colTemplate = ['varID TEXT',
        'chr TEXT',
        'pos INTEGER',
        'REF TEXT',
        'ALT TEXT',
        'QUAL TEXT',
        'FILTER TEXT',
        'INFO TEXT']
    types = [ 'TEXT' for sam in reader.samples ]
    samples = [ '\"' + sam + '\"' for sam in reader.samples  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]
    template = ','.join( colTemplate + samCol )

    ## push vcf rows into db
    with dbCon:

        ## instantiate cursor
        cur = dbCon.cursor()

        ## initialize db table
        cur.execute("DROP TABLE IF EXISTS %s" % tableName)
        cur.execute("CREATE TABLE %s(%s);" % (tableName, template))

        for rec in reader:
            ## create build a record for the variant
            varRec = [ rec.CHROM + '_' + str(rec.POS),
                rec.CHROM,
                rec.POS,
                rec.REF,
                str(rec.ALT[0]),
                rec.QUAL,
                str(rec.FILTER),
                str(rec.INFO) ]
            for sam in samples:
                varRec.append(rec.genotype(sam.strip('"'))['GT'])

            ## insert into db
            parSub = ','.join( tuple('?'*len(varRec)) )
            cur.execute("INSERT INTO %s VALUES(%s)" % (tableName, parSub), varRec)



def genoConsensus(genoList):
    '''
    Create consensus of hard call genotypes:

    @genoList: a list of genotypes in the format '0/1', with missing as '.'
    @return: return a hard call of the genotype, or '*' if no consensus

    TODO:
        TEST THIS PLEASE OH GOD TEST THIS
        Sane handling of missing data -- two missing means low qual?
    '''

    if genoList[0] == genoList[1]:
        ## two values aggree
        return genoList[0]
    elif genoList[1]==genoList[2]:
        ## two values aggree
        return genoList[1]
    elif genoList[0]==genoList[2]:
        return genoList[0]
        ## two values aggree
    else:
        ## no values aggree
        return '*'


#vcfFiles = ['data/ATLAS.merged.ontarget.chr22.SM.vcf',
#    'data/GATK.multisample.ontarget.chr22.vcf.recode.vcf',
#    'data/freebayes.chr22.multisample.ontarget.vcf.recode.vcf']


def main():

    ## parse input arguments
    parser = opt.OptionParser()
    parser.add_option('--vcf-files', dest = 'vcfFiles', action = 'store',
        nargs = 3, help = 'List of vcf files which are to be merged by consensus.')
    parser.add_option('--base-out', dest = 'baseOut', action = 'store', 
        help = 'Basename for output of .vcf and .log files.')
    (options, args) = parser.parse_args()


    ## start up database connection
    con = sql.connect('consensus.db')
    cur = con.cursor()
   

    ## push vcf files into database as tables
    for vcf in options.vcfFiles:
        fileName = os.path.basename(vcf)
        table = fileName.split('.')[0]
        db = store_vcf(vcf, table, con)



    ## find all of the tables in the db
    cur.execute("SELECT * FROM sqlite_master WHERE type='table'")
    vcfTables = [ table[1] for table in cur.fetchall() if table[1]!='consensus' ] 

    
    ## find variants and samples common to all db
    varSets = list()
    sampleSets = list()
    for table in vcfTables:
        
        cur.execute('SELECT varID FROM %s' % table )
        variants = cur.fetchall()
        
        ## add set of variants to list 
        variants = set([ x[0] for x in variants ])
        varSets.append(variants)

        ## add column names
        cur.execute("PRAGMA TABLE_INFO(%s)" % table)
        colnames = [ '\"'+col[1]+'\"' for col in cur.fetchall() ]
        samples = set(colnames[9:])
        sampleSets.append(samples) 
 
    ## set of variants common to each database        
    commonVar = reduce( lambda x,y: x.intersection(y), varSets )
 
    ## set of samples common to each database
    commonSam = reduce( lambda x,y: x.intersection(y), sampleSets )

    ## create consensus table in db
    cur.execute("DROP TABLE IF EXISTS consensus")
    colTemplate = ['varID TEXT',
        'chr TEXT',
        'pos INTEGER',
        'REF TEXT',
        'ALT TEXT',
        'QUAL TEXT',
        'FILTER TEXT',
        'INFO TEXT']
    types = [ 'TEXT' for sam in commonSam ]
    #samples = [ '\"' + sam + '\"' for sam in commonSam  ]
    samples = [ sam for sam in commonSam  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]
    template = ','.join( colTemplate + samCol )
    cur.execute("CREATE TABLE consensus(%s)" % template )
    
    ## swicth to dict cursor for fast row access
    con.row_factory = sql.Row
    cur = con.cursor()
    
    ## pull sample genotypes for each variant
    for idx, var in enumerate(commonVar):
        
        ## store genotypes for a variant for each caller
        callerGenotypes = dict()
        chr = list()
        pos = list()
        alt = list()
        ref = list()
        
        for table in vcfTables:

            cur.execute("SELECT * FROM %s WHERE varID='%s'" % (table, var))       
            row = cur.fetchone()
            callerGenotypes[table] = row
            
            ## we should check whether chr, pos, ref, alt match here
            chr.append(row['chr'])
            pos.append(row['pos'])
            ref.append(row['ALT'])
            alt.append(row['REF'])
        
        ## uniquify and check that all values that should match, match
        chr = set(chr)
        pos = set(pos)
        ref = set(ref)
        alt = set(alt)

        if len(chr)==1 and len(pos)==1 and len(ref)==1 and len(alt)==1:
            chr = chr.pop()
            pos = pos.pop()
            ref = ref.pop()
            alt = alt.pop()
        else:
            ## log these variants
            #print chr, pos, ref, alt
            #print 'variant', var, 'skipped -- more than one allele observed.'
            continue

        ## create a row of consensus genotypes for this variant
        varID = var
        qual = '-'
        filter = '-'
        info = '-'
        
        consensusRecord = [varID, chr, pos, ref, alt, qual, filter, info] 
        for sam in commonSam:
            
            ## this works but is far too slow!
            #cur.execute("SELECT %s FROM %s WHERE variant='%s'" % (sam, table, var) )
            
            ## dict access of row is much faster
            ## store genotype for each sample in consensus table
            genoSet = [ callerGenotypes[table][str(sam).strip('"')] for table in vcfTables ]
            consGeno = genoConsensus(genoSet)
            
            ## TODO :: handle missing data more intelligently
            if not consGeno:
                consGeno = './.'
            consensusRecord.append('\'' + str(consGeno) + '\'')
        
        valString = ','.join( tuple('?'*len(consensusRecord)) ) 
        cur.execute("INSERT INTO consensus VALUES(%s)" % valString, consensusRecord)
            


    ## write out a vcf file
    ## initialize vcf writer

    vcfOut = options.baseOut + '.vcf'
    vcfCon = open(vcfOut, 'w')
    ## write header
    print >> vcfCon, '##fileformat=VCFv4.0'
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for sam in commonSam: header.append(str(sam).strip('"'))
    print >> vcfCon, '\t'.join(header)

    cur.execute("SELECT * FROM consensus")
    genotypes = cur.fetchall()
    for geno in genotypes:
        ## query the record
        var = str(geno['varID'])
        chr = str(geno['chr'])
        pos = str(geno['pos'])
        ref = str(geno['REF'])
        alt = str(geno['alt'])
        sampleGeno = [ geno[str(sam).strip('"')] for sam in commonSam ]
        
        ##assemble the row
        row = [chr, pos, var, ref, alt, '.', '.', '.', '.']
        for sam in commonSam: row.append( geno[str(sam).strip('"')].strip('\'') )
        print >> vcfCon, '\t'.join(row)

if __name__ == "__main__":
    main()



