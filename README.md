galaxy.consensus
================

Consensus calling tool for cox Galaxy instance.



Introduction:
=============

This is a first pass implementation for a consensus caller that merges
genotypes from three arbitrary vcf files.


Dependencies:
=============

The major dependency that must be installed on the system is James Casbon's vcf
parser.

A cloned copy exists in lib/. Link: https://github.com/jamescasbon/PyVCF.git

You can do the standard python install by running from the PyVCF directory:

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


This should produce two files:
#### TEST.vcf
Contains consensus genotypes for variants and samples that match between all three input vcf files.
#### consensus.db
sqlite db storing vcf files for query
