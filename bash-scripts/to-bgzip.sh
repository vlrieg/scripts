#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu
module load htslib/1.3.1-gcb01
FILE=$1
bgzip ${FILE}
