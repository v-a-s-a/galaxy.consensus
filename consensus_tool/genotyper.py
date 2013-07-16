import sqlite3 as sql
import vcf as pyvcf
from consensus_functions import *

class genotyper:
  '''
  Genotyper object:

  Parses VCF files into database upon init.
  Calls consensus variants and genotypes.
  Writes consensus VCF to file.
  '''


  def __init__(self, *args, **kwargs):
    
    ## identify which callers are being used
    self.vcfTables = [ f for f in kwargs if 'Vcf' in f ]
 
    ## store VCF file locations
    self.sqliteConnection = kwargs.get('sqliteConnection', None)
    self.atlasVCF = kwargs.get('atlasVcf', None)
    self.gatkVCF = kwargs.get('gatkVcf', None)
    self.freebayesVCF = kwargs.get('freebayesVcf', None)


    ## add vcf files to database
    for branchVcf in self.vcfTables:
      self.store_vcf(vcf=kwargs.get(branchVcf), source=branchVcf)

  def constr_consensus_query(self, ncallers):
    '''
    Construct a query returning the variants concordant among *ncallers
    '''

    ## total number of tables in this run
    ntables = len(self.vcfTables)

    ## "level of consensus" -- from strict (all talbes)
    consenLevel = ntables - (ntables - ncallers)

    if ncallers > ntables:
      raise Exception('Asking for consensus among more callers than provided!')
   
    
    table1 = self.vcfTables[0]
    base = 'SELECT %s.varID FROM %s' % (table1, table1)
    extension = [ '%s ON %s.varID=%s.varID' % (table, table1, table) for table in self.vcfTables[1:] ]
    extension.insert(0, base)
    query = ' JOIN '.join(extension)

    return query

  def store_vcf(self, vcf, source):
    '''
    Inputs:
    @vcf: the file path for the vcf file to source
    @source: the name of the caller, used to identify the table
    

    Table format:
        -row is variant
        -col for varID s chr_pos
        -col for chromosome
        -col for position
        -col for REF allele
        -col for ALT allele
        -col for QUAL
        -col for sample genotype -- sample IDs are enclosed by "
    '''

    ## initialize vcf parser
    reader = pyvcf.Reader(open(vcf, 'r'))


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
    with self.sqliteConnection:

        ## instantiate cursor
        cur = self.sqliteConnection.cursor()

        ## initialize db table
        cur.execute("DROP TABLE IF EXISTS %s" % source)
        cur.execute("CREATE TABLE %s(%s);" % (source, template))

        for rec in reader:        
            
            ## skip multi-nucleotide polymorphisms
            if len(rec.ALT) > 1:
                print '\tNOT PROCESSING MULTI-ALLELIC SITES OR MNPs'
                continue

            ## create build a record for the variant
            varRec = [ rec.CHROM+'_'+str(rec.POS)+'-'+ str(rec.REF)+':'+str(rec.ALT),
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
            cur.execute("INSERT INTO %s VALUES(%s)" % (source, parSub), varRec)




  def call_consensus(self, consThresh=None):
    '''
    Reduce variant tables into consensus table.
    '''

    if consThresh is None:
      ## default consensus is among all callers input
      consThresh = len(self.vcfTables)

    cur = self.sqliteConnection.cursor()

    # get IDs of samples common to all sets
    names = list()
    for table in self.vcfTables:
        colq = cur.execute( 'SELECT * from %s' % table )
        names.append(set(map(lambda x: x[0], colq.description)))
    commonSam = reduce(lambda x,y: x.intersection(y), names )
    vcfCols = {'varID':1, 'chr':1, 'pos':1, 'REF':1, 'ALT':1, 'QUAL':1, 'FILTER':1, 'INFO':1}
    commonSam = [ x for x in commonSam if not vcfCols.get(x)  ]

    if consThresh == 3:
        ## get IDs of snps common to all sets
        varquery = 'SELECT atlas.varID FROM atlas \
                    JOIN freebayes ON atlas.varID=freebayes.varID \
                    JOIN gatk ON atlas.varID=gatk.varID'
    elif consThresh == 2:
        ## get IDs of snps in at least 2/3 sets
        varquery = 'SELECT atlas.varID FROM atlas \
                        JOIN freebayes ON atlas.varID=freebayes.varID \
                    UNION \
                    SELECT freebayes.varID FROM freebayes \
                        JOIN gatk ON freebayes.varID=gatk.varID \
                    UNION \
                    SELECT atlas.varID FROM atlas \
                        JOIN gatk ON atlas.varID=gatk.varID'

    ## pull IDs from database
    cur.execute(varquery)
    commonVar = cur.fetchall()

    print self.constr_consensus_query(3) 

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
    self.sqliteConnection.row_factory = sql.Row
    cur = self.sqliteConnection.cursor()

     ## pull sample genotypes for each variant
    for idx, var in enumerate(commonVar):
        var = var[0]

        ## store genotypes for a variant for each caller
        callerGenotypes = dict()
        infoFields = dict()

        for table in self.vcfTables:
            cur.execute("SELECT * FROM %s WHERE varID='%s'" % (table, var))
            row = cur.fetchone()
            if row:
                ## variant was observed by this algorithm
                callerGenotypes[table] = row

        
                ## grab the meta data from arbitrary caller -- all are matched on uniq var IDs
                chr = row['chr']
                pos = row['pos']
                ref = row['REF']
                alt = row['ALT']
       
                ## store QUAL scores
                infoFields[table] = row['QUAL']
            else:
                ## variant was NOT observed by the algorithm
                infoFields[table] = '-'
                callerGenotypes[table] = dict( zip(commonSam, [ None for x in commonSam]) )


        ## create a row of consensus genotypes for this variant
        varID = var
        qual = '-'
        filter = 'PASS'

        ## store QUAL scores in info matching the a specific order
        info = list()
        qualAbbrev = [ caller[0].capitalize()+'Q' for caller in self.vcfTables ]
        fieldAbbrev = dict(zip(self.vcfTables, qualAbbrev))
        for caller in fieldAbbrev.keys():
            info.append(fieldAbbrev[caller] + '=' + infoFields[caller])
        info = ';'.join(info)

        consensusRecord = [varID, chr, pos, ref, alt, qual, filter, info]
        for sam in commonSam:

            ## dict access of row is much faster
            ## store genotype for each sample in consensus table
            genoSet = [ callerGenotypes[table][str(sam).strip('"')] for table in self.vcfTables ]
            if consThresh == 3:
                genotypeField = strict_consensus(genoSet)
            elif consThresh == 2:
                genotypeField = loose_consensus(genoSet)

            ## TODO :: handle missing data more intelligently
            if not genotypeField:
                 genotypeField = './.'
            consensusRecord.append('\'' + str(genotypeField) + '\'')

        valString = ','.join( tuple('?'*len(consensusRecord)) )
        cur.execute("INSERT INTO consensus VALUES(%s)" % valString, consensusRecord)

  
  
  def make_out_vcf(commonSam):
    '''
    Only conceived of now.
    The idea is to setup the output VCF with specified parameters.
    Write records to it later, after its been initialized.
    '''
    pass

  def write_vcf(self, vcfOut):
    '''
    Write consensus VCF out to file.
    '''

    cur = self.sqliteConnection.cursor()

    names = list()
    for table in self.vcfTables:
        colq = cur.execute( 'SELECT * from %s' % table )
        names.append(set(map(lambda x: x[0], colq.description)))
    commonSam = reduce(lambda x,y: x.intersection(y), names )
    vcfCols = {'varID':1, 'chr':1, 'pos':1, 'REF':1, 'ALT':1, 'QUAL':1, 'FILTER':1, 'INFO':1}
    commonSam = [ x for x in commonSam if not vcfCols.get(x)  ]

    vcfCon = open(vcfOut, 'w')
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
       
        ## if all genotypes are missing, do not write the record out
        print set(row)
        print >> vcfCon, '\t'.join(row)




