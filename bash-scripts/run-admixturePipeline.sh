#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

module load python/3.7.4-gcb01
module load plink #pink1.9
module load vcftools
module load bcftools

infile=$1

# I updated vcf.py in the admixturePipeline dir to take a chromosome map when generating the plink files.
# Make the chrom map file here:
bcftools view ${infile}  | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt

admixturePipeline.py -m popmap.txt -v ${infile} -k 1 -K 6 -n 16 -a 0.05 -t 100 -R 100
#k smallest # of populations you want to model
#K biggest # of populations you want to model
#n number of processors for Admixture to use
#a min allele freq
#t Filter loci by thinning out any loci falling within the specified proximity to one another, measured in basepairs. 
#  (default = off, specify an integer greater than 0 to turn it on).
#R number of independant runs for each value of K