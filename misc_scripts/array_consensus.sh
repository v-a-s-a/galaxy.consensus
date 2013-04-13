#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N strict.consensus
#$ -t 1-22


time python consensus_tool/consensus.py --out ts_exomes_data/consensus.consensus.strict.chr${SGE_TASK_ID}.vcf \
    --atlas-vcf ts_exomes_data/atlas_branch/atlas_exome_chr${SGE_TASK_ID}.recode.vcf \
    --freebayes-vcf ts_exomes_data/freebayes_branch/freebayes_03-25_exome_minQ_chr${SGE_TASK_ID}.recode.vcf \
    --gatk-vcf ts_exomes_data/gatk_branch/gatk_05-05_exome_chr${SGE_TASK_ID}.recode.vcf \
    --consensus-thresh 3 \
    --db-file ~/tmp/strict.consensus.${SGE_TASK_ID}.sqlite

