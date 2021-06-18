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


#make sure there's a tabix index file
echo ">>> creating index for start file <<<"
tabix -p vcf ${NAME}.g.vcf.gz

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

#make chromosome map for generating ped file
echo ">>> making chromosome map file with bcftools <<<"
bcftools view -H ${NAME}_masked-rm.recode.chroms-only.snps.vcf  | cut -f 1 | uniq | awk '{print $0"\t"$0}' > ${NAME}.chrom-map.txt

#make ped file
echo ">>> making ped file <<<"
vcftools --vcf ${NAME}_masked-rm.recode.chroms-only.snps.vcf --plink --chrom-map ${NAME}.chrom-map.txt --out ${NAME}

#Now pick up on step 4 of the tutorial

#4. Convert to BCF and do some renaming of variants
echo ">>> rename variants with bcftools <<<"
bcftools norm -m-any --check-ref w -f ${REF} ${NAME}_masked-rm.recode.chroms-only.snps.vcf | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | bcftools norm -Ob --rm-dup both > ${NAME}.bcf
bcftools index ${NAME}.bcf

# 5 convert to PLINK format
echo ">>> convert to plink format <<<"
plink --bcf ${NAME}.bcf \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr \
      --make-bed \
      --out ${NAME}.genotypes

#6 NA (not in tutorial - artifact from previous tutorials)

#7 prune variants with MAF < 10%
echo ">>> prune variants with MAF < 10% <<<"

mkdir Pruned

plink --bfile ${NAME}.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --allow-extra-chr \
      --out Pruned/${NAME}.pruned.genotypes
      #plink2 hasn't implemented --indep flag yet

echo ">>> step 2 of pruning variants with plink <<<"      
plink --bfile ${NAME}.genotypes \
      --extract Pruned/${NAME}.pruned.genotypes.prune.in \
      --make-bed \
      --allow-extra-chr \
      --out Pruned/${NAME}.pruned.genotypes
      
#8 & 9 are merging individual chromosome files, but since I didn't split into chromosomes in the first place, I skip these steps

#10 perform pca with plink
echo ">>> performing pca with plink <<<"
plink --bfile ${NAME}.genotypes --pca --allow-extra-chr
# *** NOTE: if you don't have 50 individuals, you will get an error at this step!

#11 generate plots in R
# See tutorial https://www.biostars.org/p/335605/
