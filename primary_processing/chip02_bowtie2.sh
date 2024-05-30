#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list.txt | while read line
do
    col1="${line}"
    singularity exec bowtie2_latest.sif bowtie2 -p 32 -x Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome -1 analysis/chip/fastq/${col1}_1.paired.fq.gz -2 analysis/chip/fastq/${col1}_2.paired.fq.gz > analysis/chip/bowtie/${col1}.trimmed.sam
done


