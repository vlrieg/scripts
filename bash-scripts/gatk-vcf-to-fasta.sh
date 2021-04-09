#!/bin/bash

module load jdk/1.8.0_45-fasrc01

#PvP01 vivax reference genome
reference=/data/wraycompute/malaria/reference/vivax/PVP01.fa
ref=${reference}


#input VCF you want to convert to fasta format
vcf=$1
output=${vcf%%.vcf}.fasta

echo "vcf is ${vcf}"
echo "output is ${output}"

RUN_GATK3="java -jar /data/wraycompute/val/bin/GenomeAnalysisTK.jar"

${RUN_GATK3} -T FastaAlternateReferenceMaker -R ${ref} -o ${output} -V ${vcf}
