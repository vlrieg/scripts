#!/usr/bin/env python3

from pandas.core.dtypes import missing
from chrom_length_dict import * # /Users/vgartner/Dropbox (Duke Bio_Ea)/scripts/python-scripts/chrom_length_dict.py
import sys
import numpy as np
import pandas as pd
from operator import itemgetter
from itertools import groupby
import os.path
import csv
import collections #default dict



lmiss_file = sys.argv[1]
outfilename = lmiss_file + '_contiguous_ranges.csv'

missingness_table = pd.read_table(lmiss_file, sep = "\t")

# split SNP column 
missingness_table['POS']=missingness_table['SNP'].str.split(':').str[1]
missingness_table['REF']=missingness_table['SNP'].str.split(':').str[2]
missingness_table['ALT']=missingness_table['SNP'].str.split(':').str[3]

#drop the SNP column
missingness_table = missingness_table.drop('SNP', 1)

#rename columns
missingness_table.rename(columns={'N_GENO': 'N_TOTAL'}, inplace=True)
missingness_table.rename(columns={'F_MISS': 'FRACTION_MISSING'}, inplace=True)

#print(missingness_table.shape)
#print(missingness_table.info)
#grouped_df = missingness_table.groupby("CHR", 1)

# def missing_chunks(chr):
#     total_len = pv_chrom_length_lookup(chr)

    # number of positions with missingness per chr

    # number_missing / total_len

    # distance between missing positions ? - count?

    # return absolute number, fraction of chromosome with missing data, and some data structre with distance info?
    #return 


chrom_list = set(missingness_table['CHR'])
#print(chrom_list)

# for i in chrom_list:
#     print(i)

#print(f' overall table is: {missingness_table.info()}')


def contiguous_range(chromosome):
    chrom_df = missingness_table[missingness_table.CHR == chromosome]

    pos_list = chrom_df['POS'].tolist()
    pos_list = [int(x) for x in pos_list] # make each item in list an integer
    
    # find contiguous regions as tuples
    ranges = []
    for k,g in groupby(enumerate(pos_list),lambda x:x[0]-x[1]): #k,g == key,group
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append((group[0],group[-1]))
          
    #print(f'chromosome {chromosome} contiguous regions are: {ranges}')

    #find length of continuous regions
    results_list = []
    for tuple in ranges:
        pos_start = tuple[0]
        pos_end = tuple[1]
        range_length = pos_end - pos_start
    
        results_list.append([chromosome, pos_start, pos_end, range_length])
    return results_list


# now find the ranges for each chromosome and write to file

for each_chr in chrom_list:
    range_list_to_print = contiguous_range(each_chr)
    #print(range_list_to_print)
    """
    write to file (append)
    """

    file_exists = os.path.isfile(outfilename)
    with open(outfilename, 'a') as out_file: #'a' means results will be appended to outfilename if it exists
            headers = ['CHR', 'start', 'end', 'range_length']
            out_writer = csv.writer(out_file, delimiter=',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
            
            if not file_exists: # if true, the file doesn't exist yet -> write the header
                    out_writer.writerow(headers)

            range_len_dict = collections.defaultdict(int) #how many of each size region are there? ranges >0
            for line in range(len(range_list_to_print)):
                    chr_name = range_list_to_print[line][0]
                    pos_s = range_list_to_print[line][1]
                    pos_e = range_list_to_print[line][2]
                    rlength = range_list_to_print[line][3]
                    out_writer.writerow([chr_name, pos_s, pos_e, rlength])

                    if rlength > 0:
                        range_len_dict[rlength] += 1
            print(f'{each_chr} has the distribution of lengths:')
            for key, val in range_len_dict.items():
                print(f'{key} bp: {val} observations')            

    out_file.close()


"""
>>>>> sort the data, then:    
data = [2, 3, 4, 5, 12, 13, 14, 15, 16, 17]
ranges =[]
for k,g in groupby(enumerate(data),lambda x:x[0]-x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int,group))
    ranges.append((group[0],group[-1]))
print(ranges) -> outputs:
    [(2, 5), (12, 17)]
https://stackoverflow.com/a/2154437/10176950
"""


# for i in [unique chromosome]:
#     find the missing chunks using missing_chunks() function