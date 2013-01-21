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

from pyVCF_record import _Record, _AltRecord, _Substitution


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
        return '*/*'


#vcfFiles = ['data/ATLAS.merged.ontarget.chr22.SM.vcf',
#    'data/GATK.multisample.ontarget.chr22.vcf.recode.vcf',
#    'data/freebayes.chr22.multisample.ontarget.vcf.recode.vcf']


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
    (options, args) = parser.parse_args()


    ## start up database connection in temporary location
    con = sql.connect('/tmp/consensus.db')
    cur = con.cursor()
   

    print 'Building database of calls.'
    ## push vcf files into database as tables
    vcfFiles = {'freebayes':options.freebayesVcf,
                'gatk':options.gatkVcf,
                'atlas':options.atlasVcf}
    for vcf in vcfFiles.itervalues():
        fileName = os.path.basename(vcf)
        table = fileName.split('.')[0]
        print '\tparsing %s VCF file ...' % table
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

    print 'Calling consensus genotypes for each variant.'
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



