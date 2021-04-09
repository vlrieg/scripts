#!/bin/bash
#SBATCH --mem=60G
#SBATCH -c16
#SBATCH --job-name=bcftools-variant-call
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

# run like
# for i in $(< drc-accessions.txt ) ; do sbatch ./bcftools-variant-calling.sh ${i} ; done

module load bcftools

#REF=/data/wraycompute/malaria/reference/vivax/PVP01.fa
REF=/data/wraycompute/malaria/reference/plasmo-combined.fasta
ACC=$1

bcftools mpileup -a FORMAT/AD,FORMAT/DP -Ou ${ACC}/${ACC}.dedup.bam -f ${REF} | bcftools call -Ou    -m   --ploidy 1 | bcftools filter -Ob -i 'F_MISSING<0.1&&MAF>0.1' > ./bcftools-output/${ACC}-bcftools-called.bcf

#from https://www.biostars.org/p/9463195/
