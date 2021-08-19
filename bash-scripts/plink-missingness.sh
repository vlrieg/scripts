#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --job-name=plink-missingness

# run like ./plink-missingness.sh ethiopia-haploid-combined-joint-called
# where 'ethiopia-haploid-combined-joint-called' is the name of the vcf file WITHOUT any extensions
# though note that tabix expects the file to have the extension .g.vcf.gz

# based on this tutorial: https://www.biostars.org/p/335605/
# note this script is slightly different than my PCA script (does not removed masked regions)

module load plink
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

#extract chromosomes only (no contigs)
if [ ! -e ${NAME}.chroms-only.vcf ] && [ ! -e ${NAME}.chroms-only.snps.vcf ] ; then
  echo ">>> extract chromosomes only <<<"
  bcftools view ${NAME}.g.vcf.gz --regions-file PVP01.chroms.bed > ${NAME}.chroms-only.vcf
fi

#keep snps only
if [ ! -e ${NAME}.chroms-only.snps.vcf ] ; then
  echo ">>> keep snps only <<<"
  bcftools view -m2 -M2 -v snps ${NAME}.chroms-only.vcf > ${NAME}.chroms-only.snps.vcf 
  #https://www.biostars.org/p/141156/#141164
fi

#remove intermediate file to save storage space
if [ -e ${NAME}.chroms-only.vcf ] ; then 
  rm -rf ${NAME}.chroms-only.vcf
fi

#make chromosome map for generating ped file
if [ ! -e ${NAME}.chrom-map.txt ] ; then
  echo ">>> making chromosome map file with bcftools <<<"
  bcftools view -H ${NAME}.chroms-only.snps.vcf  | cut -f 1 | uniq | awk '{print $0"\t"$0}' > ${NAME}.chrom-map.txt
fi

#make ped file
if [ ! -e ${NAME}.ped ] || [ ! -e ${NAME}.map ] ; then
  echo ">>> making .ped and .map files <<<"
  ulimit -n 3000
  vcftools --vcf ${NAME}.chroms-only.snps.vcf --plink --chrom-map ${NAME}.chrom-map.txt --out ${NAME}
fi

# Convert to BCF and do some renaming of variants
echo ">>> rename variants with bcftools <<<"
bcftools norm -m-any --check-ref w -f ${REF} ${NAME}.chroms-only.snps.vcf | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | bcftools norm -Ob --rm-dup both > ${NAME}.bcf
bcftools index ${NAME}.bcf

# convert to PLINK format
echo ">>> convert to plink format <<<"
plink --bcf ${NAME}.bcf \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr \
      --make-bed \
      --out ${NAME}.genotypes

# find missingness before MAF pruning
echo ">>> find missingness per site and per individual PRE-PRUNING <<<"
mkdir pre-pruned
plink --bfile ${NAME}.genotypes --missing --allow-extra-chr
mv plink.lmiss pre-pruned/
mv plink.imiss pre-pruned/

echo ">>> fixing plink output format tab spacing <<<"
cat pre-pruned/plink.lmiss | tr -s ' ' '\t' > pre-pruned/plink.lmiss.tabs
cat pre-pruned/plink.imiss | tr -s ' ' '\t' > pre-pruned/plink.imiss.tabs


# prune variants with MAF < 10%
echo ">>> prune variants with MAF < 10% <<<"
mkdir Pruned

plink --bfile ${NAME}.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --allow-extra-chr \
      --out Pruned/${NAME}.pruned.genotypes

echo ">>> step 2 of pruning variants with plink <<<"      
plink --bfile ${NAME}.genotypes \
      --extract Pruned/${NAME}.pruned.genotypes.prune.in \
      --make-bed \
      --allow-extra-chr \
      --out Pruned/${NAME}.pruned.genotypes

echo ">>> find missingness per site and per individual <<<"
plink --bfile Pruned/${NAME}.pruned.genotypes --missing --allow-extra-chr

echo ">>> fixing plink output format tab spacing <<<"
cat Pruned/plink.lmiss | tr -s ' ' '\t' > Pruned/plink.lmiss.tabs
cat Pruned/plink.imiss | tr -s ' ' '\t' > Pruned/plink.imiss.tabs
