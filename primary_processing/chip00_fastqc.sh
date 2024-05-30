#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list.txt | while read line
do
    col1="${line}"
    singularity exec fastqc_0.11.9--0.sif fastqc analysis/chip/fastq/${col1}_1.fq.gz -o analysis/chip/fastqc
    singularity exec fastqc_0.11.9--0.sif fastqc analysis/chip/fastq/${col1}_2.fq.gz -o analysis/chip/fastqc
done