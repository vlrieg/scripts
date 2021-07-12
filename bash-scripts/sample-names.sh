#!/bin/bash
# run like:
# ./sample-names.sh all-samples-combined-joint-called.g.vcf > acc.txt    

# or just run command on its own (not as script) - should be very fast
# bcftools query -l all-samples-combined-joint-called.g.vcf > acc.txt

module load bcftools
INFILE=$1
bcftools query -l ${INFILE}
