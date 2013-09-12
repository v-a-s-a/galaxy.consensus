Introduction:
==============================

This is an initial implementation of a two stage voting scheme among variant calling algorithms. Given a set of VCF files produced by various algorithms, sites are selected if they are seen among all callers. Genotypes among these sites are then selected as those that match among all callers. Currently, a user can input any number of sorted VCF files, and a strict consensus of variant sites and genotypes will be generated.

Any VCF can be used as long as it can be parsed by (James Casbon's pyVCF module)[https://github.com/jamescasbon/PyVCF].

Options:
========

    usage: consensus_genotyper.py [-h] VCFS [VCFS ...]

    Find sites and genotypes which aggree among an arbitrary number of VCF files.
    
    positional arguments:
      VCFS        List of VCF files for input.

    optional arguments:
      -h, --help  show this help message and exit



Usage:
========

Test data is located in the data/ directory. The following command:

python ./consensus_tool/consensus_genotyper.py ./data/*vcf > test.output.vcf

Will take the three test files in the data directory and generate a strict consensus of sites and genotypes (i.e. 3/3 files containt the variant site, and 3/3 files agree on the genotype for a sampple at that site).

Some things to keep in mind: 
* Multi-sample VCF files are currently supported, and the output will contain only samples which are found in all input files.
* Files must be sorted by physical position. The caller works by iterating simultaneously across all input files until a matching variant record is found. If a VCF file is not sorted similarly, it is unlikely that any overlapping sites will be found.
* Missing data on the genotype level is ignored if actual genotypes are available in other VCF files. Missing data is produced only if all sites are missing, or if genotypes do not agree among all call sets.


Planned Extensions:
===================
* Operating on specific regions using the tabix index.
* Support for multi-allelic sites.
* Outputting variant sites which are discordant between callers. This is potentially interesting variation.
* The ability to specify concordance thresholds on the site and genotype level. This could be particularly helpful if one set of variants is markedly different from others, or if one is interested in finding the union of call sets rather than an intersection.
* The ability to preserve information from input VCF files. I'm thinking that it would help to specify this information in a high level configuration file. This would allow you to do things like propagate QUAL scores and compute with them downstream.

