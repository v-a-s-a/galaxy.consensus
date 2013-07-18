import sqlite3 as sql
import vcf as pyvcf
from consensus_functions import *
import itertools as iter


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
  
    ## TODO: factor out these ugly conditionals 
    base = 'SELECT %s.varID FROM %s' % (self.vcfTables[0], self.vcfTables[0])
    if consenLevel == ntables:
      ## calling the most string level of consensus
      table1 = self.vcfTables[0]
      extension = [ ' JOIN %s ON %s.varID=%s.varID' % (table, table1, table) for table in self.vcfTables[1:] ]
      extension.insert(0, base)
      query = ' '.join(extension)
    elif consenLevel == 1:
      ## allowing any variant seen anywhere
      extension = [ 'SELECT %s.varID FROM %s' % (table, table) for table in self.vcfTables[1:] ]
      extension.insert(0, base)
      query = ' UNION '.join(extension)
    elif consenLevel == 2 and ntables == 3:
      
      ## choose combinations
      combs = iter.combinations(self.vcfTables, ncallers)
      query = list()
      for combo in combs:
        ## construct select statement
        selectState = ' SELECT %s.varID FROM %s' % (combo[0], combo[0])
        ## construct join statement
        joinState = ' JOIN %s ON %s.varID=%s.varID ' % (combo[1], combo[0], combo[1])
        query.append(selectState + joinState)
      query = ' UNION '.join(query)

    else:
      raise Exception('Consensus among %i callers given %i inputs is not currently supported.' % (ncallers, ntables))

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
    ## construct genotype columns
    types = [ 'TEXT' for sam in reader.samples ]
    samples = [ '"' + sam + '"' for sam in reader.samples  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]

    ## construct DEPTH columns -- flagged as DP in most vcfs
    realTypes = [ 'REAL' for sam in reader.samples ]
    dpSamples = [ '"'+sample+'.DP'+'"' for sample in reader.samples ]
    dpCol = [ ' '.join(pair) for pair in zip(dpSamples, realTypes) ]

    ## construct PL columns
    plSamples = [ '"'+sample+'.PL'+'"' for sample in reader.samples ]
    plCol = [ ' '.join(pair) for pair in zip(plSamples, types) ]

    template = ','.join( colTemplate + samCol + dpCol + plCol )

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
                print '\tNot processing site with ALT: %s' % rec.ALT
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
            ## append genotypes
            for sam in samples:
                varRec.append(rec.genotype(sam.strip('"'))['GT'])

            ## append depths if relevant to this caller
            for sam in samples:
                try:
                    dp = getattr(rec.genotype(sam.strip('"')).data, 'DP')
                    varRec.append(dp)
                except AttributeError:
                    varRec.append('null')
            
            ## append PL values if relevant to this caller
            for sam in samples:
                try:
                    pl = str(getattr(rec.genotype(sam.strip('"')).data, 'PL'))
                    varRec.append(pl)
                except AttributeError:
                    varRec.append('null')
      
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
    ## Jesus this is ugly, and not at all extensible
    commonSam = [ x for x in commonSam if not vcfCols.get(x) and not 'PL' in x and not 'DP' in x  ]

    if consThresh == 3:
        ## get IDs of snps common to all sets
        varquery = self.constr_consensus_query(3)

    elif consThresh == 2:
        ## get IDs of snps comming to union of all pairs of sets
        varquery = self.constr_consensus_query(2)

    elif consThresh == 1:
        ## get IDs of any snp observed in the study
        varquery = self.constr_consensus_query(1)
 
    else:
        raise Exception('Consensus threshold %i is not supported!' % consThresh )

    ## pull IDs from database
    print 'CONSENSUS QUERY:\n%s' % varquery
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
    realTypes = [ 'REAL' for sam in commonSam ]
    ## genotype columns
    samples = [ '"'+sam+'"' for sam in commonSam  ]
    samCol = [ ' '.join(pair) for pair in zip(samples, types) ]

    ## average read depth columns
    adSamples = [ '"'+sample+'.AD'+'"' for sample in commonSam ]
    adCol = [ ' '.join(pair) for pair in zip(adSamples, realTypes) ]

    ## GATKs PL columns
    plSamples = [ '"'+sample+'.PL'+'"' for sample in commonSam ]
    plCol = [ ' '.join(pair) for pair in zip(plSamples, types) ]

    template = ','.join( colTemplate + samCol + adCol + plCol )
    cur.execute("CREATE TABLE consensus(%s)" % template )

    ## swicth to dict cursor for fast row access
    self.sqliteConnection.row_factory = sql.Row
    cur = self.sqliteConnection.cursor()

    ## pull sample genotypes for each variant
    for idx, var in enumerate(commonVar):
        var = var[0]

        ## store genotypes for a variant for each caller
        callerGenotypes = dict()
        callerDepths = dict()
        infoFields = dict()
        callerPLs = dict()

        for table in self.vcfTables:
            cur.execute("SELECT * FROM %s WHERE varID='%s'" % (table, var))
            row = cur.fetchone()
            if row:
                ## variant was observed by this algorithm
                callerGenotypes[table] = dict( (sam,row[sam]) for sam in commonSam )
                callerDepths[table] = dict( (sam+'.DP',row[sam+'.DP']) for sam in commonSam )
                callerPLs[table] = dict( (sam+'.PL',row[sam+'.PL']) for sam in commonSam )

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
        genotypes = list()
        adValues = list()
        plValues = list()
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
            genotypes.append('\'' + str(genotypeField) + '\'')
            #consensusRecord.append('\'' + str(genotypeField) + '\'')

            ## calculate and append the average depth metrics
            depthSet = [ callerDepths[table][sam+'.DP'] for table in self.vcfTables ]
            cleanSet = [ f for f in depthSet if f ]
            meanDP = reduce(lambda x, y: x + y, cleanSet) / len(cleanSet)
            adValues.append(meanDP)
            #consensusRecord.append(str(meanDP))

            ## append the PL metrics
            plSet = [ callerPLs[table][sam+'.PL'] for table in self.vcfTables ]
            pl = [ f.strip('[]').replace(' ','') for f in plSet if f!='none' ][0]
            plValues.append(str(pl))
            #consensusRecord.append(str(pl))

        consensusRecord = consensusRecord + genotypes + adValues + plValues
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
    commonSam = [ x for x in commonSam if not vcfCols.get(x) and not 'PL' in x and not 'DP' in x  ]

    vcfCon = open(vcfOut, 'w')
    print >> vcfCon, '##fileformat=VCFv4.0'
    print >> vcfCon, '##INFO=<ID=AQ,Number=1,Type=Float,Description="ATLAS QUAL score for variant.">'
    print >> vcfCon, '##INFO=<ID=FQ,Number=1,Type=Float,Description="Freebayes QUAL score for variant.">'
    print >> vcfCon, '##INFO=<ID=GQ,Number=1,Type=Float,Description="GATK QUAL score for variant.">'
    print >> vcfCon, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    print >> vcfCon, '##FORMAT=<ID=CN,Number=1,Type=Character,Description="Consensus status of genotype.">'
    print >> vcfCon, '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="GATK\'s Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'
    print >> vcfCon, '##FORMAT=<ID=AD,Number=1,Type=Float,Description="Average depth across callers.">'

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
        format = 'GT:CN:PL:AD'

        ## assemble the row
        row = [chr, pos, varid, ref, alt, qual, filter, info, format]
        for sam in commonSam:
            sam = str(sam).strip('"')

             ## pull genotype
            
            geno = var[sam].strip('\'')

            ## record consensus status
            if geno == '*/*':
              consensusFlag = 'F'
              geno = './.'
            else:
              consensusFlag = 'T'
            
            ## record PL
            pl = var[sam+'.PL']

            ## record average read depth
            ad = str(var[sam+'.AD'])

            ## append consensus flag
            row.append( geno+':'+consensusFlag+':'+pl+':'+ad )

       

        ## if all genotypes are missing, do not write the record out
        print >> vcfCon, '\t'.join(row)




