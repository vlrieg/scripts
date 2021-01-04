#!/bin/bash
#SBATCH --mem=60G
#SBATCH -c16
#SBATCH --job-name=joint-call
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

#run like ./joint-calling.sh population.g.vcf

#https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Discover_variants_with_GATK_-_A_GATK_Workshop_Tutorial.md
#https://github.com/kkorunes/2018_GeneConversion_Pipeline/blob/master/Scripts_SnpCallingAndRecalibration/recalibration1.2.5_GenotypeGVCFs_B.sh

module load jdk/1.8.0_45-fasrc01

# Exit immediately if any command returns a failing exit status
set -e

reference=/data/wraycompute/malaria/reference/vivax/PVP01.fa
ref=${reference}

VCF=$1 #the combined gVCF file
OUT=${VCF%%.g.vcf}-joint-called.g.vcf

echo 'VCF file is' ${VCF}
echo 'output file is' ${OUT}


GATK3_JAR=/data/wraycompute/malaria/Applications/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
RUN_GATK3="java -jar ${GATK3_JAR}"


${RUN_GATK3} -T GenotypeGVCFs -R ${ref} --variant ${VCF} -allSites -o ${OUT}
