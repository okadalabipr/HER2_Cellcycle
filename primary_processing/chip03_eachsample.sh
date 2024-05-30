#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list.txt | while read line
do
    col1="${line}"
    singularity exec samtools_1.14--hb421002_0.sif samtools view -S -h -b -F 0x4 -q 42 analysis/chip/bowtie/${col1}.trimmed.sam > analysis/chip/bowtie/${col1}.trimmed.bam
    singularity exec samtools_1.14--hb421002_0.sif samtools sort analysis/chip/bowtie/${col1}.trimmed.bam > analysis/chip/bowtie/${col1}.trimmed.sorted.bam
    singularity exec picard_2.26.6--hdfd78af_0.sif picard MarkDuplicates I=analysis/chip/bowtie/${col1}.trimmed.sorted.bam O=analysis/chip/bowtie/${col1}.uniq.proper_pairs.bam M=analysis/chip/bowtie/${col1}.uniq.proper_pairs.bam.log.txt REMOVE_DUPLICATES=true
    singularity exec samtools_1.14--hb421002_0.sif samtools sort analysis/chip/bowtie/${col1}.uniq.proper_pairs.bam -o analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bam
    singularity exec samtools_1.14--hb421002_0.sif samtools view -h analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bam > analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.sam
    singularity exec samtools_1.14--hb421002_0.sif samtools index analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bam
    singularity exec homer_4.11--pl5262h7d875b9_5.sif makeTagDirectory analysis/chip/homer/${col1} -single analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.sam
    singularity exec deeptools_3.5.1--py_0.sif bamCoverage -b analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bam -p 16 --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --binSize 1 --ignoreForNormalization chrX -o analysis/chip/bowtie/${col1}.uniq.proper_pairs.sorted.bw

done

