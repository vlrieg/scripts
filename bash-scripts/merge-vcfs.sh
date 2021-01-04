#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load vcftools

vcf-merge ../cambodia-10/cambodia-combined-joint-called.g.vcf.gz ../ethiopia-10/ethiopia-combined-joint-called.g.vcf.gz | bgzip -c > by-pop-merged.g.vcf.gz 
