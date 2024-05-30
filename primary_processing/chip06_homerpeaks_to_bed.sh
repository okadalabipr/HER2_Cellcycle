#!/bin/bash

sed '/^#/d' 12H_CHIP27AC_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 12H_CHIP27AC_histone.bed
sed '/^#/d' 16H_CTL27AC_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 16H_CTL27AC_histone.bed
sed '/^#/d' 20H_CTL27AC_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 20H_CTL27AC_histone.bed
sed '/^#/d' 24H_CTL27AC_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 24H_CTL27AC_histone.bed
sed '/^#/d' 28H_CTL27AC_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 28H_CTL27AC_histone.bed

sed '/^#/d' 12H_CHIPME1_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 12H_CHIPME1_histone.bed
sed '/^#/d' 16H_CTLME1_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 16H_CTLME1_histone.bed
sed '/^#/d' 20H_CTLME1_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 20H_CTLME1_histone.bed
sed '/^#/d' 24H_CTLME1_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 24H_CTLME1_histone.bed
sed '/^#/d' 28H_CTLME1_regions.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > 28H_CTLME1_histone.bed