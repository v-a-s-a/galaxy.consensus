class consensus_vcf:
  '''
  High level description of a VCF file.

  We want this to be able to do several things:
    1. Write a proper header given a description of fields.
    2. Correctly write records given a [branch records] and {samples, genotypes}

  '''

  @property
  def format(self):
    return self.format
  @format.setter
  def format(self, id, number, type, description):
    self.format = [ {"ID":id, "Number":number, "Type":type, "Description":'\"'+description+'\"'} ] 
  @format.getter
  def format(self):
    return self.format

  def add_format(self, id, number, type, description):
    '''
    Add data for a FORMAT field to the VCF file.
      Format is a list of dictionaries, each dict describing a format entry.
    '''
    self.format.append( {"ID":id, "Number":number, "Type":type, "Description":'\"'+description+'\"'} )


  @property
  def info(self):
    return self.info
  @info.setter
  def info(self, id, number, type, description):
    self.info = [ {"ID":id, "Number":number, "Type":type, "Description":'\"'+description+'\"'} ]
    
  def add_info(self, id, number, type, description):
    '''
    Add data for an INFO field to the VCF file.
      Info is a list of dictionaries, each dict describing an info entry.
    '''
    self.info.append( {"ID":id, "Number":number, "Type":type, "Description":'\"'+description+'\"'} )
    

  @property
  def samples(self):
    return self.samples
  @samples.setter
  def samples(self, sampleList):
    self.samples = sampleList


  def __init__(self):
    '''
    Start with a minimal description, and fill out fields from there.
    '''
  
    self.format = list()
    self.info = list()
    self.samples = list()  


  def make_var_id(self, rec):
    '''
    All variants will be uniquely identified by chr:pos:ref:alt
    '''
    return ':'.join([rec.CHROM, str(rec.POS), rec.REF, str(rec.ALT)])

  def write_header(self):
    '''
    Write a valid VCF header.
    '''
    ## FORMAT and INFO fields have a certain order to them
    order = [ 'ID', 'Number', 'Type', 'Description' ]

    ## default obligatory info
    print '##fileformat=VCFv4.1'
    ## write all INFO fields out
    for info in self.info:
      fields = ','.join([ '='.join([x, info[x]]) for x in order ]) 
      print "##INFO=<%s>" % fields
    ## write all FORMAT fields out
    for format in self.format:
      fields = ','.join([ '='.join([x, format[x]]) for x in order ])
      print "##FORMAT=<%s>" % fields
    ## write final and most important line of the header
    mainLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % '\t'.join(self.samples)
    print mainLine


  def convert_genotype(self, genotype):
    '''
    convert dosage back to alleles for writing out.
    '''
    if genotype == 0: return '0/0'
    elif genotype == 1: return '0/1'
    elif genotype == 2: return '1/1'
    else: return './.'

  def write_record(self, recordSet, genotypes):
    '''
    Write a set of genotypes/data as a row in the VCF file.
    '''
    ## fill out static fields first
    chr = recordSet[0].CHROM
    pos = str(recordSet[0].POS)
    id = self.make_var_id(recordSet[0])
    ref = recordSet[0].REF
    alt = str(recordSet[0].ALT)
    qual = '-'
    filter = 'PASS'
    info = '-'
    format = ':'.join([ x["ID"] for x in self.format ])
    
    ## put together sample fields according to FORMATs order
    genotypeFields = list()
    for sample in self.samples:
      ## TODO: include where data comes from in self.format -- this is hacky
      if genotypes[sample] == '*':
        cn = 'F'
        geno = './.'
      else:
        cn = 'T'
        geno = self.convert_genotype(genotypes[sample])
      genotypeFields.append(cn+':'+geno)
      
    ## put together the line in the VCF file
    vcfLine = [chr, pos, id, ref, alt, qual, filter, format] + genotypeFields
    print '\t'.join(vcfLine)
