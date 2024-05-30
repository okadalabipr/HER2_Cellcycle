#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

singularity exec fastqc_0.11.9--0.sif fastqc analysis/rna/fastq/*.fq.gz -o analysis/rna/fastqc

for i in "0H" "2H" "4H" "8H" "12H"
do
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_1_1.fq.gz analysis/rna/fastq/${i}_1_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_2_1.fq.gz analysis/rna/fastq/${i}_2_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_3_1.fq.gz analysis/rna/fastq/${i}_3_2.fq.gz -o trim/ramdaseq
done

singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_1_CDK4_1.fq.gz analysis/rna/fastq/16H_1_CDK4_2.fq.gz -o trim/ramdaseq
singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_1_CDK5_1.fq.gz analysis/rna/fastq/16H_1_CDK5_2.fq.gz -o trim/ramdaseq
singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_1_CDK6_1.fq.gz analysis/rna/fastq/16H_1_CDK6_2.fq.gz -o trim/ramdaseq
singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_1_DMSO_1.fq.gz analysis/rna/fastq/16H_1_DMSO_2.fq.gz -o trim/ramdaseq
singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_2_DMSO_1.fq.gz analysis/rna/fastq/16H_2_DMSO_2.fq.gz -o trim/ramdaseq
singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/16H_3_DMSO_1.fq.gz analysis/rna/fastq/16H_3_DMSO_2.fq.gz -o trim/ramdaseq

for i in "20H" "24H" "28H"
do
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_1_DMSO_1.fq.gz analysis/rna/fastq/${i}_1_DMSO_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_2_DMSO_1.fq.gz analysis/rna/fastq/${i}_2_DMSO_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_3_DMSO_1.fq.gz analysis/rna/fastq/${i}_3_DMSO_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_1_CDK4_1.fq.gz analysis/rna/fastq/${i}_1_CDK4_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_2_CDK4_1.fq.gz analysis/rna/fastq/${i}_2_CDK4_2.fq.gz -o trim/ramdaseq
    singularity exec trim_galore_latest.sif trim_galore --trim1 --paired analysis/rna/fastq/${i}_3_CDK4_1.fq.gz analysis/rna/fastq/${i}_3_CDK4_2.fq.gz -o trim/ramdaseq
done