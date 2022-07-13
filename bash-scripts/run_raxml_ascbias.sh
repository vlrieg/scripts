#!/usr/bin/env bash
#SBATCH --job-name=raxml
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

infile=$1
raxml-ng --all --msa ${infile} --tree pars{10} --bs-tree 100 --model GTR+ASC_LEWIS --prefix ${infile}.LEWIS --threads 16
