#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --job-name=missingness

# run like ./gtovcf_then_find-missingness.sh ethiopia-haploid-combined-joint-called
# where 'ethiopia-haploid-combined-joint-called' is the name of the vcf file WITHOUT any extensions
# though note that tabix expects the file to have the extension .g.vcf.gz

# based on this tutorial: https://www.biostars.org/p/335605/
# note this script is slightly different than my PCA script (does not removed masked regions)

module load tabix
module load vcftools/0.1.13-gcb01
module load htslib/1.3.1-gcb01
module load bcftools/1.3.1-gcb01

NAME=$1
REF=/data/wraycompute/malaria/reference/vivax/PVP01.fa
set -e

#create a tabix index file if not already present
if [ ! -e ${NAME}.g.vcf.gz.tbi ] ; then
  echo ">>> creating index for start file <<<"
  tabix -p vcf ${NAME}.g.vcf.gz
fi

#remove masked regions
if [ ! -e ${NAME}_masked-rm.recode.vcf ] && [ ! -e ${NAME}_masked-rm.recode.vcf.gz ] ; then 
  echo ">>> removing masked regions <<<"
  vcftools --gzvcf ${NAME}.g.vcf.gz --out ${NAME}_masked-rm --exclude-bed PVP01.genome.mask.sorted.bed --recode --keep-INFO-all
  #out file will look like
  #ethiopia-haploid-combined-joint-called_masked-rm.recode.vcf
fi

if [ ! -e ${NAME}_masked-rm.recode.vcf.gz ] ; then
  bgzip ${NAME}_masked-rm.recode.vcf
  bcftools index ${NAME}_masked-rm.recode.vcf.gz
fi

#extract chromosomes only (no contigs)
if [ ! -e ${NAME}_masked-rm.recode.chroms-only.vcf ] && [ ! -e ${NAME}_masked-rm.recode.chroms-only.vcf.gz ] && [ ! -e ${NAME}_masked-rm.recode.chroms-only.snps.vcf ] ; then
  echo ">>> extract chromosomes only <<<"
  bcftools view ${NAME}_masked-rm.recode.vcf.gz --regions-file PVP01.chroms.bed > ${NAME}_masked-rm.recode.chroms-only.vcf
fi

#keep snps only
if [ ! -e ${NAME}_masked-rm.recode.chroms-only.snps.vcf ] && [ ! -e ${NAME}_masked-rm.recode.chroms-only.snps.vcf.gz ]; then
  echo ">>> keep snps only <<<"
  bcftools view -m2 -M2 -v snps ${NAME}_masked-rm.recode.chroms-only.vcf > ${NAME}_masked-rm.recode.chroms-only.snps.vcf 
  #https://www.biostars.org/p/141156/#141164
fi

#remove intermediate file to save storage space
if [ -e ${NAME}_masked-rm.recode.chroms-only.vcf ] ; then 
  rm -rf ${NAME}_masked-rm.recode.chroms-only.vcf
fi

# vcftools missingness
if [ ! -e ${NAME}_masked-rm.recode.chroms-only.snps.vcf.gz ] ; then
  bgzip ${NAME}_masked-rm.recode.chroms-only.snps.vcf
  tabix -p vcf ${NAME}_masked-rm.recode.chroms-only.snps.vcf.gz
fi

# individual missingness
vcftools --gzvcf ${NAME}_masked-rm.recode.chroms-only.snps.vcf.gz --missing-indv
# site missingness
vcftools --gzvcf ${NAME}_masked-rm.recode.chroms-only.snps.vcf.gz --missing-site
