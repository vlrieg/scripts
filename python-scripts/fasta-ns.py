#!/usr/bin/env python3
#run this file with two arguments - the fasta file you want counted and the acceptable percentage of Ns:
# ./fasta-ns.py cambodia-api-seqs.fasta 30

import sys
import argparse
from Bio import SeqIO #to parse the fasta file
from collections import defaultdict

parser = argparse.ArgumentParser(description="Extract sequences which do not have too many Ns")
parser.add_argument("file",help="FASTA file containing sequences to filter")
parser.add_argument("percent",type=float,default=30,help="Maximum percentags of Ns allowed")
args = parser.parse_args()

fasta_file_name = args.file

fasta_file=None
if fasta_file_name == '-':
    fasta_file=sys.stdin
else:
    fasta_file=open(fasta_file_name,"r")

fasta_list = list(SeqIO.parse(fasta_file, "fasta"))

percentage_acceptable_Ns = args.percent


def calc_n_percentage(fastas, percent):
    #fastas is a list of seq objects
    #percent is the limit for percentage of Ns per sequence
    # (e.g. for fewer than 30% Ns in a sequence allwed, this argument = 30)
    
    fasta_len = len(fastas) #the number of items in the list
    
    output_dict = defaultdict(lambda: 'Missing Sequence')
    reject_count = 0
    
    for seq in fastas:
        seq_len = len(seq)
        seq_id = seq.id
        
        n_count = 0
        gap_count = 0
        letter_count = 0
        
        s = str(seq.seq)
        n_count = s.count('N') + s.count('n') + s.count('?')
        gap_count = s.count('-')
        letter_count = seq_len - n_count - gap_count
        
        #calculate percent Ns
        percent_n = (n_count + gap_count )/ seq_len * 100
        
        #keep only seqs with fewer than some percentage of Ns
        if percent_n < percent:
            output_dict[seq_id] = seq
            print(f"KEEP   {seq_id}: {percent_n:0.1f}% Ns", file=sys.stderr)
        else:
            print(f"Reject {seq_id}: {percent_n:0.1f}% Ns", file=sys.stderr)
            reject_count += 1
    
    results = [output_dict, reject_count, fasta_len]
    return(results)



# Call the function & store in variable
output = calc_n_percentage(fasta_list, percentage_acceptable_Ns)

#store each result in its own variable
output_dict = output[0]
too_many_ns = output[1]
total_no_ids = output[2]

okay_ns = total_no_ids - too_many_ns #how many sequences are kept

#write desired fasta sequences to new output file
for key, val in output_dict.items():
    SeqIO.write(val, sys.stdout, "fasta")

#print output message
print('----------------------------------------------------------------------------------------------------------------------------', file=sys.stderr)
print(f'{okay_ns} of {total_no_ids} sequences have enough data, while {too_many_ns} sequences had too many Ns to include in analysis.' , file=sys.stderr)
print(f'{okay_ns} fasta sequences printed to new file.', file=sys.stderr)
print('----------------------------------------------------------------------------------------------------------------------------', file=sys.stderr)
