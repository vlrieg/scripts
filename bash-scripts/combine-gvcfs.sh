#!/bin/bash

# This script is supposed to be based on https://software.broadinstitute.org/gatk/documentation/article.php?id=6483
# Run like `sbatch ./combine-gvcfs.sh southamerica`

module load bwa
module load samtools
module load jdk/1.8.0_45-fasrc01

# Exit immediately if any command returns a failing exit status
set -e

#PvP01 vivax reference genome
reference=/data/wraycompute/vdp5/reference_data/PVP01.fasta
ref=${reference}

prefix=$1

RUN_GATK3="java -jar /home/vdp5/src/GenomeAnalysisTK.jar"

function CombineGVCFs ()
{
${RUN_GATK3} -T CombineGVCFs -R ${ref} --variant ${prefix}_input.list -o ${prefix}-cohort.g.vcf 
}
CombineGVCFs #call the function

