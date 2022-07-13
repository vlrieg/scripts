#!/usr/bin/env bash
#SBATCH --job-name=raxml
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

infile=$1
felsenstein_file="${infile}.felsenstein"

#read in felsenstein value
fels_val=$(<${felsenstein_file})

raxml-ng --all --msa ${infile} --tree pars{10} --bs-tree 100 --model GTR+G+ASC_FELS{${fels_val}} --prefix ${infile}.FELS --threads 16
