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
from consensus_functions import *




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
        help = 'Location of file for sqlite db')
    (options, args) = parser.parse_args()


    ## start up database connection in temporary location
    con = sql.connect(options.dbFile)
    cur = con.cursor()
   

    print 'Building database of calls.'
    ## push vcf files into database as tables
    vcfFiles = {'freebayes':options.freebayesVcf,
                'gatk':options.gatkVcf,
                'atlas':options.atlasVcf}
    vcfTables = ['freebayes', 'gatk', 'atlas']
    for table in vcfTables:
        vcf = vcfFiles[table]
        print '\tparsing %s VCF file ...' % table
        db = store_vcf(vcf, table, con)


    ## get IDs of samples common to all sets
    names = list()
    for table in vcfTables:
        colq = cur.execute( 'SELECT * from %s' % table )
        names.append(set(map(lambda x: x[0], colq.description)))
    commonSam = reduce(lambda x,y: x.intersection(y), names )
    vcfCols = {'varID':1, 'chr':1, 'pos':1, 'REF':1, 'ALT':1, 'QUAL':1, 'FILTER':1, 'INFO':1}
    commonSam = [ x for x in commonSam if not vcfCols.get(x)  ]

    ## get IDs of snps common to all sets
    varquery = 'SELECT a.varID FROM\
            atlas AS a\
            JOIN freebayes AS f ON a.varID=f.varID\
            JOIN gatk as g ON a.varID=g.varID'
    cur.execute(varquery)
    commonVar = cur.fetchall() 

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
    samples = [ '"' + sam + '"' for sam in commonSam  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]
    template = ','.join( colTemplate + samCol )
    cur.execute("CREATE TABLE consensus(%s)" % template )
    
    ## swicth to dict cursor for fast row access
    con.row_factory = sql.Row
    cur = con.cursor()


    
    ## pull sample genotypes for each variant
    for idx, var in enumerate(commonVar):
        var = var[0]
        
        ## store genotypes for a variant for each caller
        callerGenotypes = dict()
        chr = list()
        pos = list()
        alt = list()
        ref = list()
        infoFields = dict()
        
        for table in vcfTables:
            cur.execute("SELECT * FROM %s WHERE varID='%s'" % (table, var))       
            row = cur.fetchone()
            callerGenotypes[table] = row
            
            ## we should check whether chr, pos, ref, alt match here
            chr.append(row['chr'])
            pos.append(row['pos'])
            ref.append(row['ALT'])
            alt.append(row['REF'])

            ## store QUAL scores
            infoFields[table] = row['QUAL']
        
        ## TODO this should be guaranteed upstream!
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
        filter = 'PASS'

        ## store QUAL scores in info matching the a specific order
        info = list()
        fieldAbbrev = {'atlas':'AQ', 'freebayes':'FQ', 'gatk':'GQ'}
        for caller in fieldAbbrev.keys():
            info.append(fieldAbbrev[caller] + '=' + infoFields[caller])
        info = ';'.join(info)
        
        consensusRecord = [varID, chr, pos, ref, alt, qual, filter, info] 
        for sam in commonSam:
            
            ## dict access of row is much faster
            ## store genotype for each sample in consensus table
            genoSet = [ callerGenotypes[table][str(sam).strip('"')] for table in vcfTables ]
            genotypeField = genoConsensus(genoSet)
            
            ## TODO :: handle missing data more intelligently
            if not genotypeField:
                 genotypeField = './.'
            consensusRecord.append('\'' + str(genotypeField) + '\'')
        
        valString = ','.join( tuple('?'*len(consensusRecord)) ) 
        cur.execute("INSERT INTO consensus VALUES(%s)" % valString, consensusRecord)
            


    ## write out a vcf file
    ## initialize vcf writer

    vcfOut = options.baseOut
    vcfCon = open(vcfOut, 'w')
    ## write header
    print >> vcfCon, '##fileformat=VCFv4.0'
    print >> vcfCon, '##INFO=<ID=AQ,Number=1,Type=Float,Description="ATLAS QUAL score for variant.">'
    print >> vcfCon, '##INFO=<ID=FQ,Number=1,Type=Float,Description="Freebayes QUAL score for variant.">'
    print >> vcfCon, '##INFO=<ID=GQ,Number=1,Type=Float,Description="GATK QUAL score for variant.">'
    print >> vcfCon, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    print >> vcfCon, '##FORMAT=<ID=CN,Number=1,Type=Character,Description="Consensus status of genotype.">'
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for sam in commonSam: header.append(str(sam).strip('"'))
    print >> vcfCon, '\t'.join(header)

    cur.execute("SELECT * FROM consensus")
    consensusVar = cur.fetchall()
    for var in consensusVar:
        ## query variant record
        varid = str(var['varID'])
        chr = str(var['chr'])
        pos = str(var['pos'])
        ref = str(var['REF'])
        alt = str(var['alt'])
        qual = '.'
        filter = '.'
        info = str(var['INFO']) 
        format = 'GT:CN'
        #sampleGeno = [ var[str(sam).strip('"')] for sam in commonSam ]
        
        ##assemble the row
        row = [chr, pos, varid, ref, alt, qual, filter, info, format]
        for sam in commonSam:
            ## pull genotype
            geno = var[str(sam).strip('"')].strip('\'')

            ## record consensus status
            if geno == '*/*':
              consensusFlag = 'F'
              geno = './.'
            else:
              consensusFlag = 'T'
 
            row.append( geno + ':' + consensusFlag )
        print >> vcfCon, '\t'.join(row)

if __name__ == "__main__":
    main()



