#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load bcftools

bcftools view -s SANRU_9_HC_ATTACTCG-CCTATCCT_S5_L001 2021-09-07-global-combined-jointcalled.g.vcf.gz > post-jc_SANRU_9.g.vcf

