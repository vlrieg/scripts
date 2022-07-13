#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

source /gpfs/fs1/home/vag15/miniconda3/bin/activate /gpfs/fs1/home/vag15/miniconda3/envs/vivax
module load tabix

vcf=$1
popfile=$2

pixy --stats pi dxy fst --vcf ${vcf} --populations ${popfile} --window_size 100 --output_prefix global_${popfile%%.txt}_organelle_pixy --chromosomes 'LT635626,LT635627'
