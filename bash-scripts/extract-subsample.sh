#!/usr/bin/env bash
#SBATCH --job-name=bcft_view
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load bcftools
module load htslib/1.3.1-gcb01
module load tabix

VCF=$1
SUBSAMPLE=$2

bcftools view -S ${SUBSAMPLE} ${VCF} > ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}.vcf

bgzip ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}.vcf
tabix -p vcf ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}.vcf.gz
