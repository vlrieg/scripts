#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --job-name=vcftools-missingness


module load vcftools/0.1.13-gcb01


infile=$1

# individual missingness
vcftools --gzvcf ${infile} --missing-indv
mv out.imiss ${infile}-out.imiss


# site missingness
vcftools --gzvcf ${infile} --missing-site
mv out.lmiss ${i}-out.lmiss
