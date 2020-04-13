#!/usr/bin/env python3

import sys
import csv
from Bio import SeqIO #to parse the fasta file
from collections import defaultdict
 
location_info_file = sys.argv[1]

fasta_file = sys.argv[2]
fasta_list = list(SeqIO.parse(fasta_file, "fasta"))
fasta_len = len(fasta_list)

#Dictionary for storing SRR### (key) and location (value)
location_dict = defaultdict(lambda: 'Location Missing')

with open(location_info_file, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',', quotechar='|')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            #This is just the column headers
            line_count +=1
        else:
            location_dict[row[0]] = row[1]



#Dictionary for storing SRR###-LT63562# (key) and sequence (value)
fasta_dict = defaultdict(lambda: 'Something Missing')

for i in range(fasta_len):
    seq = fasta_list[i].seq
    seq_id = fasta_list[i].id

    fasta_dict[seq_id] = seq


#compare between dictionaries & write to new file
output_dict = defaultdict(lambda: 'output missing')

file = open(f"{fasta_file}-LocationUpdated.fasta", "w")
key_count = 0
key1_count = 0

for key, val in fasta_dict.items():
    key_count += 1
    for key1, val1 in location_dict.items():
        key1_count += 1
        if key1 in key: 
            file.write(f">{val1}-{key}\n{val}\n")
        else:
            pass

file.close()