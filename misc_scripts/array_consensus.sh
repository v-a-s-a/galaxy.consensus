#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N arr.consensus
#$ -t 1-22


time python consensus.py --out ts_exomes_data/consensus.chr${SGE_TASK_ID}.vcf \
    --atlas-vcf ts_exomes_data/atlas_exome_chr${SGE_TASK_ID}.recode.vcf \
    --freebayes-vcf ts_exomes_data/freebayes_03-25_exome_minQ_chr${SGE_TASK_ID}.recode.vcf \
    --gatk-vcf ts_exomes_data/gatk_05-05_exome_chr${SGE_TASK_ID}.recode.vcf \
    --db-file ~/tmp/consensus.${SGE_TASK_ID}.sqlite

