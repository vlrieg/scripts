#!/usr/bin/env python3
#to subsample singletons file for testing remove-masked-regions-singletonfile.py

import random
import sys
import csv

singfile = sys.argv[1]
header_list = ['CHROM', 'POS', 'SINGLETON/DOUBLETON', 'ALLELE', 'INDV']


with open(singfile, mode='r') as singleton_file:
     singleton_table = csv.reader(singleton_file, delimiter = '\t', quotechar='|')

     for sing_obs in singleton_table: #for each singleton observation
        if sing_obs == header_list: #don't waste time looping thru header a bunch
            print(*sing_obs, sep="\t")
        else:
            x = random.uniform(0,1)

            if x < 0.05:
                print(*sing_obs, sep="\t")
            else:
                pass