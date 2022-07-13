#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --job-name=indv-vcf

module load bcftools
module load vcftools
module load tabix
module load htslib/1.3.1-gcb01

for sample in `bcftools query -l chr-renamed_SORTED-2021-09-07-global-combined-jointcalled_masked-rm.recode.chroms-only.snps.20percmiss_site-then-indv.NO_UNKNOWNS_OR_POLYCLONAL_pub-pv-only.pruned.genotypes-recode01-pvp01chroms.vcf` ; do vcf-subset --exclude-ref -c $sample chr-renamed_SORTED-2021-09-07-global-combined-jointcalled_masked-rm.recode.chroms-only.snps.20percmiss_site-then-indv.NO_UNKNOWNS_OR_POLYCLONAL_pub-pv-only.pruned.genotypes-recode01-pvp01chroms.vcf  > filtered-for-admix-${sample}.vcf ; done
