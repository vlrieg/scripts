#!/usr/bin/env python3

from Bio import SeqIO #to parse the fasta file
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Calculate the number of segregating sites",
                                 epilog= "Example: Cambodia-apicoplast.fasta --reference PATH/references/PVP01.fasta --chromosome LT635626")


parser.add_argument("file", help="File with locus sequences in fasta format")
parser.add_argument("--reference", help="The reference genome in fasta format")
parser.add_argument("--chromosome", help="The chromosome from which to extract the consensus sequence")
args = parser.parse_args()

#reference
def load_reference_seqs(filename):
    reference = {}
    for chr_record in SeqIO.parse(filename,"fasta"):
        reference[chr_record.name] = chr_record
    return reference

reference_seqs = load_reference_seqs(args.reference)
reference = reference_seqs[args.chromosome] #reference sequence for your chromosome/locus of interest
# print(f"The reference sequence is {len(reference)} bp long")

# for i in range(len(reference)):
#     print(f"position {i} is {reference[i]}")


#data into a list by sample
locus_fasta = args.file
locus_list = list(SeqIO.parse(locus_fasta, "fasta"))
locus_len = len(locus_list) #the number of samples in the dictionary
#print(f"The the number of samples is {locus_len}")

# how to get the sequence on its own
# for i in range(locus_len):
#     print(locus_list[i].seq) #prints the sequence for that sample


for i in range(locus_len):
    #print(f"The length of sequence {i} is {len(locus_list[i].seq)}")
    assert len(locus_list[i].seq) == len(reference), "Sample sequence is not the same length as the reference"

column_alt = defaultdict(int)
# compare your sample sequences to the reference
for sample in range(locus_len):
    num_match = defaultdict(int)
    num_alt = defaultdict(int)
    num_n = defaultdict(int)
    for base in range(len(reference)):
        #print(f"{sample}: reference: {reference[base]} sample: {locus_list[sample].seq[base]}")
        
        if locus_list[sample].seq[base] != 'N' and locus_list[sample].seq[base] != '-':
            if reference[base] == locus_list[sample].seq[base]:
                #print(reference[base])
                num_match[base] = 1
            else:
                #print(f"ref: {reference[base]} alt: {locus_list[sample].seq[base]}")
                num_alt[base] = 1
                column_alt[base] = column_alt[base] + 1
        else:
            num_n[base] = 1

    counted = len(num_match) + len(num_alt) + len(num_n)
    total = len(reference)
    assert counted == total, "Number of positions counted does not match the number of positions in the reference!"
    #print(f"{counted} bases were counted out of {total} bases in the locus for sample {sample}")
    #print(f"{len(num_alt)} bases didn't match reference for sample {sample}")


singletons = 0
non_singletons = 0
for position, count in column_alt.items():
    #print(f"There are {count} alt alleles at position {position}.")
    if count == 1:
        singletons+=1
    else:
        non_singletons+=1

assert singletons + non_singletons == len(column_alt), "SNP counting has gone wrong!"
print(f"There are {singletons} singletons in this population.\nThere are {non_singletons} SNPs that are found in more than one sample.")



#### command used for testing ####
#./segregating-sites.py Cambodia-api-Ns-calculated.fasta --reference ~/Dropbox\ \(Duke\ Bio_Ea\)/references/PVP01.fasta --chromosome LT635626
#### script originally in directory: ####
# /Users/vgartner/Dropbox\ (Duke\ Bio_Ea)/coalescent/mtDNA/PythonScripts/SegregatingSites/segregating-sites.py
