#!/usr/bin/env python3
# This script takes as input a text file with more than one fasta file
# and outputs one file per accession number
# note this output will go to the directory where the script is saved, not where input data is

from Bio import SeqIO 
import sys

fasta_file_name = sys.argv[1]

fasta_file=None
if fasta_file_name == '-':
    fasta_file=sys.stdin
else:
    fasta_file=open(fasta_file_name,"r")

fasta_list = list(SeqIO.parse(fasta_file, "fasta"))

for seq in fasta_list:
    seq_id = seq.id
    sequence = seq.seq
    
    with open(f'{seq_id}.fasta', 'w') as out_file:
        out_file.write(f'>{seq_id}\n')
        out_file.write(f'{sequence}')
    
    out_file.close()

