#!/usr/bin/env python3

# This script compares two lists of SRA accession numbers ('original list' - SRA accession numbers already downloaded 
# and 'new list' with a selection of SRA accession numbers you might be interested in downloading if not dups).
# This script also uses the corresponding SRA run table downloaded from the Run Selector for the samples of interest - HOWEVER
# to make things more simple, I just copied the  "Run" and "geo_loc_name" columns and put in new .csv file manually. This .csv
# file is the third argument passed to the script.
# This script prints the unique (not yet downloaded) SRA run accession numbers and their corresponding geo-locations to a new file.

import sys
import csv

file1 = sys.argv[1] #original list
file2 = sys.argv[2] #new list
runtable = sys.argv[3] #includes geoloc data


# Read the two text files into lists
with open(file1, 'rb') as original_file:
    lines1 = original_file.read().decode('utf-8').splitlines()

with open(file2, 'rb') as new_file:
    lines2 = new_file.read().decode('utf-8').splitlines()


# Remove duplicates from the lists before comparing to each other
lines1 = list(dict.fromkeys(lines1))
lines2 = list(dict.fromkeys(lines2))


# Compare original and new list of SRA acession numbers

unique_accession = []

for line in lines2:
    if line not in lines1:
        unique_accession.append(line)

#duplicates = set(lines1).intersection(lines2)



# Print out unique accession numbers to new file

# text_file = open("Accession-output.txt", "w")
# text_file.write('\n'.join(unique_accession))
# text_file.close()


# Look up these unique accession numbers in the run table
# to extract geoloc info
accession_geoloc = open("accession_numbers_with_geoloc.txt", "w")

with open(runtable, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',', quotechar='|')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print_this = ", ".join(row)
            #accession_geoloc.write(f'{print_this} \n')
            line_count +=1
        else:
            for accession in unique_accession:
                if accession == row[0]:
                    # select only the "Run" and "geo_loc_name" columns and put in new .csv file manually
                    accession_geoloc.write(f'{row[0]}\t{row[1]}\n')
                    line_count += 1
                else: 
                    pass


accession_geoloc.close()