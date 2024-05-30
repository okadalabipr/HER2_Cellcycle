#!/bin/bash
#PBS -q MEDIUM
#PBS -l select=1:ncpus=8

for i in "0H" "2H" "4H" "8H" "12H"
do

    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_1_1_val_1.fq.gz trim/ramdaseq/${i}_1_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_1_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_2_1_val_1.fq.gz trim/ramdaseq/${i}_2_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_2_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_3_1_val_1.fq.gz trim/ramdaseq/${i}_3_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_3_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1

done

singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_1_CDK4_1_val_1.fq.gz trim/ramdaseq/16H_1_CDK4_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_1_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_1_CDK5_1_val_1.fq.gz trim/ramdaseq/16H_1_CDK5_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_2_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_1_CDK6_1_val_1.fq.gz trim/ramdaseq/16H_1_CDK6_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_3_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1

singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_1_DMSO_1_val_1.fq.gz trim/ramdaseq/16H_1_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_1_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_2_DMSO_1_val_1.fq.gz trim/ramdaseq/16H_2_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_2_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/16H_3_DMSO_1_val_1.fq.gz trim/ramdaseq/16H_3_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/16H_3_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1

for i in "20H" "24H" "28H"
do
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_1_DMSO_1_val_1.fq.gz trim/ramdaseq/${i}_1_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_1_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_2_DMSO_1_val_1.fq.gz trim/ramdaseq/${i}_2_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_2_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_3_DMSO_1_val_1.fq.gz trim/ramdaseq/${i}_3_DMSO_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_3_DMSO_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_1_CDK4_1_val_1.fq.gz trim/ramdaseq/${i}_1_CDK4_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_1_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_2_CDK4_1_val_1.fq.gz trim/ramdaseq/${i}_2_CDK4_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_2_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
    singularity exec star_2.7.9a--h9ee0642_0.sif STAR --genomeDir STAR_index/GRCh38 --runThreadN 4 --readFilesIn trim/ramdaseq/${i}_3_CDK4_1_val_1.fq.gz trim/ramdaseq/${i}_3_CDK4_2_val_2.fq.gz --outFileNamePrefix analysis/rna/star/${i}_3_CDK4_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
done



singularity exec subread_latest.sif featureCounts -p -B -t transcripts -g gene_id -a Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf -o analysis/rna/premrnatest.txt analysis/rna/star/*.bam
singularity exec subread_latest.sif featureCounts -p -B -t transcripts -g gene_id -a Hgencode.v40.chr_patch_hapl_scaff.annotation.gtf -o analysis/rna/premrnatest.txt analysis/rna/star/*.bam