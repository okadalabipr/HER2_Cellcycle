#!/bin/bash
#PBS -q SMALL
#PBS -l select=1:ncpus=1

source .bash_profile

cd analysis/chip/homer/merge/histone

findMotifsGenome.pl H3K27hh16HCTLINH_enhancer_me1_corr_sg.bed hg38 H3K27hh16HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh16HCTLINH_enhancer_me1_corr_sl.bed hg38 H3K27hh16HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh20HCTLINH_enhancer_me1_corr_sg.bed hg38 H3K27hh20HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh20HCTLINH_enhancer_me1_corr_sl.bed hg38 H3K27hh20HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh24HCTLINH_enhancer_me1_corr_sg.bed hg38 H3K27hh24HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh24HCTLINH_enhancer_me1_corr_sl.bed hg38 H3K27hh24HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh28HCTLINH_enhancer_me1_corr_sg.bed hg38 H3K27hh28HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh28HCTLINH_enhancer_me1_corr_sl.bed hg38 H3K27hh28HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1



findMotifsGenome.pl H3K27Ac_homerhistone_16HCTLINH_enhancer_me1.bed hg38 H3K27Ac_homerhistone_16HCTLINH_enhancer_me1_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_20HCTLINH_enhancer_me1.bed hg38 H3K27Ac_homerhistone_20HCTLINH_enhancer_me1_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_24HCTLINH_enhancer_me1.bed hg38 H3K27Ac_homerhistone_24HCTLINH_enhancer_me1_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_28HCTLINH_enhancer_me1.bed hg38 H3K27Ac_homerhistone_28HCTLINH_enhancer_me1_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1


findMotifsGenome.pl H3K27hh16HCTLINH_promoter_nme1ts250_sg.bed hg38 H3K27hh16HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh16HCTLINH_promoter_nme1ts250_sl.bed hg38 H3K27hh16HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh20HCTLINH_promoter_nme1ts250_sg.bed hg38 H3K27hh20HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh20HCTLINH_promoter_nme1ts250_sl.bed hg38 H3K27hh20HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh24HCTLINH_promoter_nme1ts250_sg.bed hg38 H3K27hh24HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh24HCTLINH_promoter_nme1ts250_sl.bed hg38 H3K27hh24HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27hh28HCTLINH_promoter_nme1ts250_sg.bed hg38 H3K27hh28HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27hh28HCTLINH_promoter_nme1ts250_sl.bed hg38 H3K27hh28HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1

findMotifsGenome.pl H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250.bed hg38 H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250.bed hg38 H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250.bed hg38 H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
findMotifsGenome.pl H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250.bed hg38 H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250_fasta_maskhomer/ -dumpFasta -noknown -nomotif -mask -size given -p 1
