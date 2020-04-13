#!/usr/bin/env Rscript

library(seqinr)
setwd("/mnt/hardac/gpfs/fs1/data/wraycompute/val/beast2/taylor-escalante-reproduction")
args = commandArgs(trailingOnly=TRUE) #for script version only
filename <- args[1]
sample_size <- args[2]

#Read in fasta file
dnafile <- read.fasta(file = "pvivax_mtDNA_locationFormat.muscle.fasta", as.string = TRUE)

#Get accession names
names <- getName(dnafile)

#Select random subset of samples by name
#### # Total number of samples in the fasta file
#### #total_num = length(names)
#### # Select how many you want in the subsample
#### #subset_num = total_num / 2

# Generate random subsample list
rand_subsample <- sample(names, size=sample_size, replace=F)


#Write new fasta file containing only the randomly generated subset of samples
# New object with only the subset of samples picked above
seqs_towrite <- dnafile[c(which(names(dnafile) %in% rand_subsample))]

# Make sure the working directory is correct
setwd("/mnt/hardac/gpfs/fs1/data/wraycompute/val/beast2/taylor-escalante-reproduction")

# Write to file
outfile_name <- paste("subsample_of_", sample_size, ".fasta", sep="")
filename <- write.fasta(sequences = seqs_towrite, names=names(seqs_towrite), file.out = outfile_name, open="w")