#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

cat analysis/chip/chip_list_for_homer_diff.txt | while read line
do
    col1=`echo ${line} | cut -d , -f 1`
    col2=`echo ${line} | cut -d , -f 2`
    singularity exec homer_4.11--pl5262h7d875b9_5.sif findPeaks analysis/chip/homer/merge/${col1} -style histone -L 0 -fdr 0.00005 -o analysis/chip/homer/merge/${col1}_${col2}_regions.txt -i analysis/chip/homer/merge/${col2}
    singularity exec homer_4.11--pl5262h7d875b9_5.sif findPeaks analysis/chip/homer/merge/${col2} -style histone -L 0 -fdr 0.00005 -o analysis/chip/homer/merge/${col2}_${col1}_regions.txt -i analysis/chip/homer/merge/${col1}
done