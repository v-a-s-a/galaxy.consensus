
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
        return '*/*'
