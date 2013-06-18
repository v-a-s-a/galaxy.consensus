
import vcf as vcfMod


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
        -col for varID is chr_pos_ref_alt
        -col for chromosome
        -col for position
        -col for REF allele
        -col for ALT allele
        -col for QUAL
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
    samples = [ '"' + sam + '"' for sam in reader.samples  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]
    template = ','.join( colTemplate + samCol )

    ## push vcf rows into db
    with dbCon:

        ## instantiate cursor
        cur = dbCon.cursor()

        ## initialize db table
        cur.execute("DROP TABLE IF EXISTS %s" % tableName)
        cur.execute("CREATE TABLE %s(%s);" % (tableName, template))

       # for rec in reader:        

        ## while loop for debugging record iter in reader
        while True:
        
    
            try:
                rec = reader.next()
            except StopIteration:
                break
            except ValueError:
                tableSize = cur.execute("SELECT varID from %s" % tableName)
                tableSize = str(len(tableSize.fetchall()))
                print '# of variants inserted into %s: %s' % (tableName, tableSize)
                
            
            ## skip multi-nucleotide polymorphisms
            if len(rec.ALT) > 1:
                continue

 
            ## create build a record for the variant
            varRec = [ '_'.join([rec.CHROM, str(rec.POS), rec.REF, str(rec.ALT[0])]),
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
                

def is_het(g):
    '''
    Check if hardcall genotype is heterozygous
    '''
    if g=='1/0' or g=='0/1': return True
    else: return False

def is_hom_ref(g):
    '''
    Check if hardcall genotype is homozygous reference
    '''
    if g=='0/0': return True
    else: return False

def is_hom_alt(g):
    '''
    Check if hardcall genotype if homozygous alternative
    '''
    if g=='1/1': return True
    else: return False

def loose_consensus(genoList):
    '''
    Call consensus if genotype is concordant among 2/3 callers.

    @genoList: a list of genotypes in the format '0/1', with missing as '.'
    @return: return a hard call of the genotype. */* if no consensus.

    '''

    if genoList[0] == genoList[1]:
        ## two values aggree
        return genoList[0]
    elif genoList[1] == genoList[2]:
        ## two values aggree
        return genoList[1]
    elif genoList[0] == genoList[2]:
        return genoList[0]
        ## two values aggree
    else:
        ## no values aggree
        return '*/*'


def strict_consensus(genoList):
    '''
    Call consensus if genotype is concordant among 3/3 callers.

    @genoList: a list of genotypes in the format '0/1', with missing as '.'
    @return: return a hard call of the genotype. */* if no consensus.

    '''

    if reduce(lambda x,y: x and y, map(is_het, genoList)):
        ## all genotypes are heterozygous
        return '0/1'
    elif reduce(lambda x,y: x and y, map(is_hom_ref, genoList)):
        ## all genotypes are homozygous reference
        return '0/0'
    elif reduce(lambda x,y: x and y, map(is_hom_alt, genoList)):
        ## all genotypes are homozygous alternative
        return '1/1'        
    else:
        ## no genotypes are concordant
        return '*/*'


