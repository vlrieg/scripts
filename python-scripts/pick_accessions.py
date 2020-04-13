#!/usr/bin/env python3

import sys
import random

file = sys.argv[1]

with open(file, 'rb') as accession_file:
    accession_list = accession_file.read().decode('utf-8').splitlines()

# ten_percent = int(len(accession_list) * .1) #718 is 10% of 7187 for Sanger_Pv_Accessions.txt
# five_percent = int(len(accession_list) * .05) #359 is ~5% of 7187 for Sanger_Pv_Accessions.txt
# print(f'{five_percent} is ~5% of {len(accession_list)}')

random_list = random.choices(population=accession_list, k=50) # Note: random.choices samples WITH replacement
print('\n'.join(random_list))
