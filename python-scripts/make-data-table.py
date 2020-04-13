#!/usr/bin/env python3

#Input data was downloaded from https://www.ncbi.nlm.nih.gov/nuccore by:
#search for accessions of interest on https://www.ncbi.nlm.nih.gov/sites/batchentrez
#click "Retrieve records for # UID(s)" link
#click Send To > click Complete Record > click File > Format = Summary > create file

#filenames
input_file = '/Users/vgartner/Dropbox (Duke Bio_Ea)/Miscellaneous/sars-cov-2/nuccore_result.txt'
output_file = '/Users/vgartner/Dropbox (Duke Bio_Ea)/Miscellaneous/sars-cov-2/covid-19_metadata.csv'

#Step 1. Parse data file by sample
with open(input_file, "r") as infile:
    data = infile.read().splitlines() #there are trailing spaces but will deal with these later

data_index = 0
results = []

while data_index < len(data):
    sample_line_index = 0
    sample = list(range(3))

    while sample_line_index < 3:
        index = data_index + sample_line_index
        sample[sample_line_index] = data[index]
        #data_index is the "true" line of the file
        #sample_line_index is the line within the sample
        sample_line_index += 1
    
    results.append(sample)
    data_index += 4 #each sample is a group of 4 lines


#print(results[0])
#['1. Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome', '29,903 bp linear RNA ', 'MN908947.3 GI:1798172431']


#Step 2. Reformat each line to keep only the following information and store in new list
#results[0]: location
#results[1]: sequence length (remove comma)
#results[2]: accession beginning with 'M'

results_reformatted = []

import re #grep

#iterate through each sample
for sample in range(len(results)):
    sample_reformatted = list(range(3))
    ##### location #####
    sample_name = results[sample][0]

    #remove leading info
    leading_len = len('. Severe acute respiratory syndrome coronavirus 2 isolate ')
    start_val = sample_name.index('. ') + leading_len #keep everything after first part of title
    sample_name = sample_name[start_val:]

    #remove trailing info
    end_val = sample_name.index(',')
    sample_name = sample_name[:end_val]

    ##### sequence length #####
    sample_length = results[sample][1]

    #remove trailing info
    len_end_val = sample_length.index(' bp')
    sample_length = sample_length[:len_end_val]

    #convert from string to int
    sample_length = int(sample_length.replace(',', ''))

    ##### accession #####
    accession = results[sample][2]
    gi_index = accession.index(' GI') #delete the GI accession at the end of the line
    accession = accession[:gi_index]

    sample_reformatted = [sample_name, sample_length, accession]
    results_reformatted.append(sample_reformatted)

infile.close()

#Step 3. Write to CSV file
import csv

with open(output_file, 'w') as out_file:
    out_writer = csv.writer(out_file, delimiter=',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)

    #write header line:
    out_writer.writerow(['accession', 'info', 'sample length (bp)'])

    for sample in range(len(results_reformatted)):
        ref_accession = results_reformatted[sample][2]
        ref_sample_name = results_reformatted[sample][0]
        ref_sample_length = results_reformatted[sample][1]
        out_writer.writerow([ref_accession, ref_sample_name, ref_sample_length])


out_file.close()