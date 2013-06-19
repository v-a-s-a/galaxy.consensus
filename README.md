Introduction:
==============================

This is an initial implementation of a two stage voting scheme among variant calling algorithms. Given a set of VCF files produced by various algorithms, sites are selected if they are seen among all callers. Genotypes among these sites are then selected as those that match among all callers.


Supported Callers:

  - [GATK Unified Genotyper] (http://www.broadinstitute.org/gatk/)
  - [Freebayes](https://github.com/ekg/freebayes)
  - [Atlas-SNP2](http://sourceforge.net/p/atlas2/wiki/Atlas-SNP/)




    Options:
      -h, --help            show this help message and exit
      --out                 File path base for output of .vcf and .log files.
      --atlas-vcf           Location of ATLAS vcf file for consensus.
      --gatk-vcf            Location of GATK vcf file for consensus.
      --freebayes-vcf       Location of freebayes vcf file for consensus.
      --db-file             Location of file for sqlite db.




Usage:
========

Test data is located in the data/ directory.

To run the test data, try:

    consensus.py --out TEST \
      --atlas-vcf data/atlas.small.vcf \
      --freebayes-vcf data/freebayes.small.vcf \
      --gatk-vcf data/gatk.small.vcf \
      --db-file mktemp


This produces a file: TEST.vcf. This Contains consensus genotypes for variants and samples that match between all
three input vcf files.

