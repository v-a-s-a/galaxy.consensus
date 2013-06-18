#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N loose.exon.consensus
#$ -t 1-22


time python consensus_tool/consensus.py \
    --out ts_exomes_data/consensus.loose.exonBed.atlas1.4.3.allsamples.flat.chr${SGE_TASK_ID}.vcf \
    --atlas-vcf ts_exomes_data/atlas_branch/atlas_exome_bed_v1.4.3_allsamples_flat_chr${SGE_TASK_ID}.vcf.recode.vcf \
    --freebayes-vcf ts_exomes_data/freebayes_branch/freebayes_03-25_exome_bed_chr${SGE_TASK_ID}.recode.vcf \
    --gatk-vcf ts_exomes_data/gatk_branch/gatk_05-05_exome_bed_multi_chr${SGE_TASK_ID}.vcf.recode.vcf \
    --consensus-thresh 2 \
    --db-file ~/tmp/loose.exon.consensus.${SGE_TASK_ID}.sqlite

