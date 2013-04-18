'''
This is going to be messy.

String together vcftools, plink, and python to produce a set of QC metrics on a
vcf.
'''

import subprocess as sp
import gzip as gz
import optparse as opt

def smart_open(f):
    '''
    Open compressed files if compressed
    '''
    if f.endswith('.gz'):
        return gz.open(f)
    else:
        return open(f)

def smart_vcftools(call):
    '''
    Analog to smart open for calling vcftools with a gzipped vcf
    '''
    for idx, element in enumerate(call):
        if element.endswith('.vcf.gz'):
            print element
            call[idx-1] = '--gzvcf'
    return call

def main():

    ## parse command line arguments
    parser = opt.OptionParser()
    parser.add_option('--vcf', dest = 'vcfFile', action = 'store', 
        help = 'File path for vcf for which to generate QC metrics.')
    parser.add_option('--intermed-base', dest = 'intermedBase', action = 'store', 
        help = 'File path base for storage of intermediate files.')
    (options, args) = parser.parse_args()

    vcfFile = options.vcfFile
    intermedBase = options.intermedBase

    devnull = open('/dev/null', 'w')
    sampleKeyFile = '/nas40t2/tsExoms/Sample Info/Scharf_master_sequencing_sample_key.csv'
    pedFile = '/nas40t2/tsExoms/Sample Info/updated_Scharf_Pedigree.csv'



    ## identify the total number of variants in the VCF file
    zgrepCall = ['zgrep', '-v', '\'#\'', vcfFile]
    wcCall = ['wc', '-l']
    zgrep = sp.Popen( zgrepCall, stdout = sp.PIPE )
    wc = sp.Popen( wcCall, stdin = zgrep.stdout, stdout = sp.PIPE )
    zgrep.stdout.close()
    nsnp = wc.communicate()[0]


    ## get the ts/tv ratio
    tstvCall = ['vcftools', '--vcf', vcfFile, '--TsTv', nsnp, '--out',
        intermedBase + '.tstv']
    tstvCall = smart_vcftools(tstvCall)
    #print ' '.join(tstvCall)
    #sp.call(tstvCall)



    ## convert vcf to plink formatted files
    convertCall = ['vcftools', '--vcf', vcfFile, '--plink', '--recode', '--out',
        intermedBase]
    convertCall = smart_vcftools(convertCall)
    #sp.call(convertCall)
    
    ## extract qual scores from INFO field
    qualCall = ['vcftools', '--vcf', vcfFile,
    '--get-INFO', 'AQ', '--get-INFO', 'GQ', '--get-INFO', 'FQ',
    '--out', intermedBase]
    sp.call(smart_vcftools(qualCall))

    ## store mapping of SM to RG
    rg_sm = dict()
    sm_rg = dict()
    fconn = smart_open(sampleKeyFile)
    header = fconn.next()
    for line in fconn:
       line = line.split(',')
       rg_sm[line[1]] = line[0]
       sm_rg[line[0]] = line[1]  
    fconn.close()

    ## add pedigree information to plink ped file file
    mothers = dict()
    fathers = dict()
    families = dict()
    sexes = dict()
    with smart_open(pedFile) as fconn:
        header = enumerate(fconn.next().split(','))
        index = dict( ( (x[1], x[0]) for x in header ) )
        for line in fconn:
            line = line.split(',')
            if rg_sm.get(line[index['Mother']]):
                mother = rg_sm[ line[index['Mother']] ]
            else:
                mother = line[index['Mother']]

            if rg_sm.get(line[index['Father']]):
                father = rg_sm[ line[index['Father']]]
            else:
                father = line[index['Father']]
            
            if mother=='Not Sequenced':
                mother = '0'
            if father=='Not Sequenced':
                father = '0'

            fathers[line[index['SM_tag (Local_ID)']]] = father
            mothers[line[index['SM_tag (Local_ID)']]] = mother

            if rg_sm.get(line[index['Individ']]):
                family = line[index['Family']]
            else:
                family = line[index['SM_tag (Local_ID)']] 

            families[line[index['SM_tag (Local_ID)']]] = family

            sexes[line[index['SM_tag (Local_ID)']]] = line[index['Sex']]


    newPlinkFile = intermedBase + '.pedigree.ped'
    newPlinkConn = open(newPlinkFile, 'w')
    with smart_open(intermedBase + '.ped') as plinkFile:
        for line in plinkFile:
            line = line.strip().split('\t')
            ## convert ID from RG to SM if necessary
            if rg_sm.get(line[1]):
              iid = rg_sm[line[1]]
              line[1] = iid
            else:
              iid = line[1]
            if fathers.get(iid):
                line[2] = fathers[iid]
            if mothers.get(iid):
                line[3] = mothers[iid]
            if families.get(iid):
                line[0] = families[iid]
            if sexes.get(iid):
                line[4] = sexes[iid] 
            ## add a fake phenotype
            line[5] = '2'
            print >> newPlinkConn, '\t'.join(line)
    newPlinkConn.close()

    ## calculate mendelian errors
    mendelCall = ['plink',
        '--ped', newPlinkFile,
        '--map', intermedBase + '.map',
        '--mendel',
        '--allow-no-sex',
        '--out', intermedBase]
    #sp.call(mendelCall)# stdout = devnull, stderr = devnull)

    ## generate minor allele frequency distribution
    mafCall = ['plink',
    '--ped', newPlinkFile,
    '--map', intermedBase + '.map',
    '--freq',
    '--allow-no-sex',
    '--out', intermedBase]
    #sp.call(mafCall)

    
    ## get the missingness rates
    missCall = ['plink',
    '--ped', newPlinkFile,
    '--map', intermedBase + '.map',
    '--missing',
    '--out', intermedBase]
    #sp.call(missCall)


    ## get hardy weinberg estimates
    hweCall = ['plink',
    '--ped', newPlinkFile,
    '--map', intermedBase + '.map',
    '--hardy',
    '--out', intermedBase]
    #sp.call(hweCall)


    devnull.close()



if __name__ == '__main__':
    main()
