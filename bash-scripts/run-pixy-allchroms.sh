#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

source /gpfs/fs1/home/vag15/miniconda3/bin/activate /gpfs/fs1/home/vag15/miniconda3/envs/vivax
module load tabix

vcf=$1
popfile=$2

pixy --stats pi dxy fst --vcf ${vcf} --populations ${popfile} --window_size 50000 --chromosomes 'LT635612,LT635613,LT635614,LT635615,LT635616,LT635617,LT635618,LT635619,LT635620,LT635621,LT635622,LT635623,LT635624,LT635625' --output_prefix global_${popfile%%.txt}_allchroms_pixy
