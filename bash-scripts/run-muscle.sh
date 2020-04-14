#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

#run on hardac with:
#sbatch run-muscle.sh sequences.fasta

module load muscle/3.8.31-fasrc01 

file=$1

muscle -in ${file} -out ${file%%.fasta}-muscle.fasta
