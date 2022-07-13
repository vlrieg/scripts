#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

# run like ./generate-plink-pca.sh ethiopia-haploid-combined-joint-called
# where 'ethiopia-haploid-combined-joint-called' is the name of the vcf file WITHOUT any extensions
# though note that tabix expects the file to have the extension .g.vcf.gz

# based on this tutorial: https://www.biostars.org/p/335605/
# note: PCA wants at least 50 individuals to impute allele frequencies from!

module load plink
module load tabix
module load vcftools
module load bcftools

NAME=$1
REF=/data/wraycompute/malaria/reference/vivax/PVP01.fa
set -e


#remove masked regions
echo ">>> removing masked regions <<<"
vcftools --gzvcf ${NAME}.g.vcf.gz --out ${NAME}_masked-rm --exclude-bed PVP01.genome.mask.sorted.bed --recode --keep-INFO-all
#out file will look like
#ethiopia-haploid-combined-joint-called_masked-rm.recode.vcf

#extract chromosomes only (no contigs)
echo ">>> extract chromosomes only <<<"
bgzip ${NAME}_masked-rm.recode.vcf
bcftools index ${NAME}_masked-rm.recode.vcf.gz
bcftools view ${NAME}_masked-rm.recode.vcf.gz --regions-file PVP01.chroms.bed > ${NAME}_masked-rm.recode.chroms-only.vcf

#keep snps only
echo ">>> keep snps only <<<"
bcftools view -m2 -M2 -v snps ${NAME}_masked-rm.recode.chroms-only.vcf > ${NAME}_masked-rm.recode.chroms-only.snps.vcf 
#https://www.biostars.org/p/141156/#141164