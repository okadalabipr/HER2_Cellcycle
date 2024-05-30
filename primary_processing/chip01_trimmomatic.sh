#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list.txt | while read line
do
    col1="${line}"
    singularity exec trimmomatic_v0.39.sif trimmomatic PE -threads 16 -phred33 analysis/chip/fastq/${col1}_1.fq.gz analysis/chip/fastq/${col1}_2.fq.gz analysis/chip/fastq/${col1}_1.paired.fq.gz analysis/chip/fastq/${col1}_1.unpaired.fq.gz analysis/chip/fastq/${col1}_2.paired.fq.gz analysis/chip/fastq/${col1}_2.unpaired.fq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:2
done