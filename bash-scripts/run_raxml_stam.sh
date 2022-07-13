#!/usr/bin/env bash
#SBATCH --job-name=raxml
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

infile=$1
stam_file="${infile}.stamatakis"

# read in A/C/G/T values
stam_vals=$(<${stam_file})

# subset for one value per variable
stam_a=$(echo $stam_vals | cut -d' ' -f 1)
stam_c=$(echo $stam_vals | cut -d' ' -f 2)
stam_g=$(echo $stam_vals | cut -d' ' -f 3)
stam_t=$(echo $stam_vals | cut -d' ' -f 4)


raxml-ng --all --msa ${infile} --tree pars{10} --bs-tree 100 --model GTR+G+ASC_STAM{${stam_a}/${stam_c}/${stam_g}/${stam_t}} --prefix ${infile}.STAM --threads 16
