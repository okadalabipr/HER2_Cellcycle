#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list_for_merge.txt | while read line
do
    col1=`echo ${line} | cut -d , -f 1`
    col2=`echo ${line} | cut -d , -f 2`
    col3=`echo ${line} | cut -d , -f 3`
    singularity exec samtools_1.14--hb421002_0.sif samtools merge analysis/chip/bowtie/merge/${col3}.merged.bam analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bam analysis/chip/bowtie/${col2}.uniq.proper_pairs.sorted.bam
    singularity exec samtools_1.14--hb421002_0.sif samtools sort analysis/chip/bowtie/merge/${col3}.merged.bam -o analysis/chip/bowtie/merge/${col3}.merged.sorted.bam
    singularity exec samtools_1.14--hb421002_0.sif samtools view -h analysis/chip/bowtie/merge/${col3}.merged.sorted.bam > analysis/chip/bowtie/merge/${col3}.merged.sorted.sam
    singularity exec samtools_1.14--hb421002_0.sif samtools index analysis/chip/bowtie/merge/${col3}.merged.sorted.bam
    singularity exec homer_4.11--pl5262h7d875b9_5.sif makeTagDirectory analysis/chip/homer/merge/${col3} -single analysis/chip/bowtie/merge/${col3}.merged.sorted.sam
    singularity exec deeptools_3.5.1--py_0.sif bamCoverage -b analysis/chip/bowtie/merge/${col3}.merged.sorted.bam -p 16 --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --binSize 1 --ignoreForNormalization chrX -o analysis/chip/bowtie/merge/${col3}.merged.sorted.bw
done