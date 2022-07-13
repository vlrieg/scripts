#!/usr/bin/env bash
#SBATCH --job-name=iqtree
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

infile=$1
iqtree -s ${infile} -m GTR+ASC -nt AUTO
