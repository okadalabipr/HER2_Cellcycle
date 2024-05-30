#!/bin/bash



bedtools intersect -v -a 12H_CHIP27AC_histone.bed -b tss250.bed > 12H_CHIP27AC_histone_ntss.bed
bedtools intersect -wa -a 12H_CHIP27AC_histone_ntss.bed -b 12H_CHIPME1_histone.bed > 12H_CHIP27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 16H_CTL27AC_histone.bed -b tss250.bed > 16H_CTL27AC_histone_ntss.bed
bedtools intersect -wa -a 16H_CTL27AC_histone_ntss.bed -b 16H_CTLME1_histone.bed > 16H_CTL27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 20H_CTL27AC_histone.bed -b tss250.bed > 20H_CTL27AC_histone_ntss.bed
bedtools intersect -wa -a 20H_CTL27AC_histone_ntss.bed -b 20H_CTLME1_histone.bed > 20H_CTL27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 24H_CTL27AC_histone.bed -b tss250.bed > 24H_CTL27AC_histone_ntss.bed
bedtools intersect -wa -a 24H_CTL27AC_histone_ntss.bed -b 24H_CTLME1_histone.bed > 24H_CTL27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 28H_CTL27AC_histone.bed -b tss250.bed > 28H_CTL27AC_histone_ntss.bed
bedtools intersect -wa -a 28H_CTL27AC_histone_ntss.bed -b 28H_CTLME1_histone.bed > 28H_CTL27AC_histone_enhancer_me1.bed

bedtools intersect -wa -a 16H_INH27AC_histone.bed -b tss250.bed > 16H_INH27AC_histone_tss250.bed
bedtools intersect -v -a 16H_INH27AC_histone_tss250.bed -b 16H_INHME1_histone.bed > 16H_INH27AC_histone_promoter_nme1ts250.bed
bedtools intersect -wa -a 20H_INH27AC_histone.bed -b tss250.bed > 20H_INH27AC_histone_tss250.bed
bedtools intersect -v -a 20H_INH27AC_histone_tss250.bed -b 20H_INHME1_histone.bed > 20H_INH27AC_histone_promoter_nme1ts250.bed
bedtools intersect -wa -a 24H_INH27AC_histone.bed -b tss250.bed > 24H_INH27AC_histone_tss250.bed
bedtools intersect -v -a 24H_INH27AC_histone_tss250.bed -b 24H_INHME1_histone.bed > 24H_INH27AC_histone_promoter_nme1ts250.bed
bedtools intersect -wa -a 28H_INH27AC_histone.bed -b tss250.bed > 28H_INH27AC_histone_tss250.bed
bedtools intersect -v -a 28H_INH27AC_histone_tss250.bed -b 28H_INHME1_histone.bed > 28H_INH27AC_histone_promoter_nme1ts250.bed

bedtools intersect -v -a 16H_INH27AC_histone.bed -b tss250.bed > 16H_INH27AC_histone_ntss.bed
bedtools intersect -wa -a 16H_INH27AC_histone_ntss.bed -b 16H_INHME1_histone.bed > 16H_INH27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 20H_INH27AC_histone.bed -b tss250.bed > 20H_INH27AC_histone_ntss.bed
bedtools intersect -wa -a 20H_INH27AC_histone_ntss.bed -b 20H_INHME1_histone.bed > 20H_INH27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 24H_INH27AC_histone.bed -b tss250.bed > 24H_INH27AC_histone_ntss.bed
bedtools intersect -wa -a 24H_INH27AC_histone_ntss.bed -b 24H_INHME1_histone.bed > 24H_INH27AC_histone_enhancer_me1.bed
bedtools intersect -v -a 28H_INH27AC_histone.bed -b tss250.bed > 28H_INH27AC_histone_ntss.bed
bedtools intersect -wa -a 28H_INH27AC_histone_ntss.bed -b 28H_INHME1_histone.bed > 28H_INH27AC_histone_enhancer_me1.bed


cat 16H_CTL27AC_histone_promoter_nme1ts250.bed 16H_INH27AC_histone_promoter_nme1ts250.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250.bed
rm -rf tmp1.bed tmp2.bed

cat 20H_CTL27AC_histone_promoter_nme1ts250.bed 20H_INH27AC_histone_promoter_nme1ts250.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250.bed
rm -rf tmp1.bed tmp2.bed

cat 24H_CTL27AC_histone_promoter_nme1ts250.bed 24H_INH27AC_histone_promoter_nme1ts250.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250.bed
rm -rf tmp1.bed tmp2.bed

cat 28H_CTL27AC_histone_promoter_nme1ts250.bed 28H_INH27AC_histone_promoter_nme1ts250.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250.bed
rm -rf tmp1.bed tmp2.bed


annotatePeaks.pl H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/16H_CTL27AC ../chip/homer_dir/16H_INH27AC > H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/20H_CTL27AC ../chip/homer_dir/20H_INH27AC > H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/24H_CTL27AC ../chip/homer_dir/24H_INH27AC > H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/28H_CTL27AC ../chip/homer_dir/28H_INH27AC > H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250_signal.txt



cat 16H_CTL27AC_histone_enhancer_me1.bed 16H_INH27AC_histone_enhancer_me1.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_16HCTLINH_enhancer_me1.bed
rm -rf tmp1.bed tmp2.bed

cat 20H_CTL27AC_histone_enhancer_me1.bed 20H_INH27AC_histone_enhancer_me1.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_20HCTLINH_enhancer_me1.bed
rm -rf tmp1.bed tmp2.bed

cat 24H_CTL27AC_histone_enhancer_me1.bed 24H_INH27AC_histone_enhancer_me1.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_24HCTLINH_enhancer_me1.bed
rm -rf tmp1.bed tmp2.bed

cat 28H_CTL27AC_histone_enhancer_me1.bed 28H_INH27AC_histone_enhancer_me1.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > H3K27Ac_homerhistone_28HCTLINH_enhancer_me1.bed
rm -rf tmp1.bed tmp2.bed


annotatePeaks.pl H3K27Ac_homerhistone_16HCTLINH_enhancer_me1.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/16H_CTL27AC ../chip/homer_dir/16H_INH27AC > H3K27Ac_homerhistone_16HCTLINH_enhancer_me1_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_20HCTLINH_enhancer_me1.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/20H_CTL27AC ../chip/homer_dir/20H_INH27AC > H3K27Ac_homerhistone_20HCTLINH_enhancer_me1_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_24HCTLINH_enhancer_me1.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/24H_CTL27AC ../chip/homer_dir/24H_INH27AC > H3K27Ac_homerhistone_24HCTLINH_enhancer_me1_signal.txt

annotatePeaks.pl H3K27Ac_homerhistone_28HCTLINH_enhancer_me1.bed hg38 -norm 1000000 -size given -d ../chip/homer_dir/28H_CTL27AC ../chip/homer_dir/28H_INH27AC > H3K27Ac_homerhistone_28HCTLINH_enhancer_me1_signal.txt