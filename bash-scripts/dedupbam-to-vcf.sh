#!/bin/bash
#SBATCH --mem=60G
#SBATCH --job-name=vcf
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load java/1.8.0_45-fasrc01
ACC=$1

java -jar /data/wraycompute/malaria/Applications/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/wraycompute/malaria/reference/plasmo-combined.fasta -I ${ACC}/${ACC}.dedup.bam -ploidy 1 -o ${ACC}/${ACC}-unfiltered-haploid.vcf 
