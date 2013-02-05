galaxy.consensus Introduction:
==============================

This is a first pass implementation for a consensus caller that merges
output VCF files from ATLAS SNP2, freebayes, and GATK Unified Genotyper.

This repo contains a tool appropriate for incorporation into Galaxy.

Options:
========
    -h, --help                        show this help message and exit
    --out=BASEOUT                     File path base for output of .vcf and .log files.
    --atlas-vcf=ATLASVCF              Location of ATLAS vcf file for consensus.
    --gatk-vcf=GATKVCF                Location of GATK vcf file for consensus.
    --freebayes-vcf=FREEBAYESVCF      Location of freebayes vcf file for
consensus.


Dependencies:
=============

The major dependency that must be installed on the system is James Casbon's vcf
parser.

A cloned copy exists in lib/, otherwise:
https://github.com/jamescasbon/PyVCF.git

You can do the standard python install by running from ./lib/PyVCF:

    python setup.py build
    python setup.py install


Example:
========

Test data is located in the data/ directory.

To run the test data, try:

    consensus.py --base-out TEST \
      --vcf-files data/ATLAS.merged.ontarget.chr22.SM.vcf \
      data/freebayes.chr22.multisample.ontarget.vcf.recode.vcf \
      data/GATK.multisample.ontarget.chr22.vcf.recode.vcf


This should produce a file:
#### TEST.vcf
Contains consensus genotypes for variants and samples that match between all
three input vcf files.

