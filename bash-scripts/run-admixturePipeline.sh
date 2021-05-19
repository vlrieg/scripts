#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

module load python/3.7.4-gcb01
module load plink #pink1.9
module load vcftools

infile=$1

admixturePipeline.py -m popmap.txt -v ${infile} -k 1 -K 6 -n 16 -a 0.05 -t 100 -R 50
