Introduction:
==============================

An initial implementation of a consensus calling algorithm, and a variant quality/summary report.


Currently, the tools supports the following callers:
  - GATK Unified Genotyper
  - Freebayes
  - Atlas-SNP2




Dependencies:
=============

The major dependency that must be installed on the system is James Casbon's vcf
parser.

A cloned copy exists in lib/, otherwise:
https://github.com/jamescasbon/PyVCF.git



Options:
========
    -h, --help                        show this help message and exit
    --out=BASEOUT                     File path base for output of .vcf and .log files.
    --atlas-vcf=ATLASVCF              Location of ATLAS vcf file for consensus.
    --gatk-vcf=GATKVCF                Location of GATK vcf file for consensus.
    --freebayes-vcf=FREEBAYESVCF      Location of freebayes vcf file for consensus.
    --db-file=DBFILE                  Location of file for sqlite db


Example:
========

Test data is located in the data/ directory.

To run the test data, try:

    consensus.py --out TEST \
      --atlas-vcf data/ATLAS.merged.ontarget.chr22.SM.vcf \
      --freebayes-vcf data/freebayes.chr22.multisample.ontarget.vcf.recode.vcf \
      --gatk-vcf data/GATK.multisample.ontarget.chr22.vcf.recode.vcf \
      --db-file mktemp


This produces a file:
#### TEST.vcf
Contains consensus genotypes for variants and samples that match between all
three input vcf files.

