#!/bin/bash
#PBS -q SMALL
#PBS -l select=1:ncpus=1

source .bash_profile

cd ichikawa-sensei/chip/homer/merge/histone


####promoter

sea --oc 16HCTLINHsg_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh16HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 16HCTLINHsl_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh16HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 20HCTLINHsg_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh20HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 20HCTLINHsl_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh20HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 24HCTLINHsg_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh24HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 24HCTLINHsl_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh24HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 28HCTLINHsg_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh28HCTLINH_promoter_nme1ts250_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 28HCTLINHsl_motif_meme_homerhistone_promoter_nme1ts250_maskhomer --p H3K27hh28HCTLINH_promoter_nme1ts250_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme



####enhancer

sea --oc 16HCTLINHsg_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh16HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_16HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 16HCTLINHsl_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh16HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_16HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 20HCTLINHsg_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh20HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_20HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 20HCTLINHsl_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh20HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_20HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 24HCTLINHsg_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh24HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_24HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 24HCTLINHsl_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh24HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_24HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme

sea --oc 28HCTLINHsg_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh28HCTLINH_enhancer_me1_corr_sg_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_28HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
sea --oc 28HCTLINHsl_motif_meme_homerhistone_enhancer_me1_corr_maskhomer --p H3K27hh28HCTLINH_enhancer_me1_corr_sl_fasta_maskhomer/target.fa --n H3K27Ac_homerhistone_28HCTLINH_enhancer_me1_fasta_maskhomer/target.fa --m ../../../../../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme