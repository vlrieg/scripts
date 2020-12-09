#!/usr/bin/env python3

# run like:
# ./remove-masked-regions-singletonfile.py ../../masks/mask_regions.csv by-pop-chroms-snps-only.singletons

# the two tables look like:

# 1. mask_regions.csv:
# chr	exclude_start	exclude_end
# LT635612	0	116540
# LT635612	903590	1021663
# LT635613	0	100045

# 2. singleton file:
# CHROM	POS	SINGLETON/DOUBLETON	ALLELE	INDV
# LT635615	105	S	A	SRR2316105  
# LT635615	127	S	C	SRR2316105
# LT635615	141	S	G	SRR2316875
#   note - this file was generated using vcftools:
#   $ vcftools --gzvcf file.g.vcf.gz --singletons


import sys
import csv
from collections import defaultdict

maskfile = sys.argv[1]
singfile = sys.argv[2]

masks_for_chrom = {}
with open(maskfile, mode='r') as mask_file:
    mask_table = csv.DictReader(mask_file, delimiter = ',')
    for row in mask_table:
        chrom = row['chr']
        start = int(row['exclude_start'])
        end   = int(row['exclude_end'])
        if chrom not in masks_for_chrom:
            masks_for_chrom[chrom] = []
        masks_for_chrom[chrom].append([start,end])


keep_dict = defaultdict(str)
keep_count = 0
header_list = ['CHROM', 'POS', 'SINGLETON/DOUBLETON', 'ALLELE', 'INDV']


# Read singleton file into a list and then compare to mask table
# only keep lines that don't fall within masked regions
# also keep all lines that are from apicoplast and mitochondria
with open(singfile, mode='r') as singleton_file:
    singleton_table = csv.reader(singleton_file, delimiter = '\t', quotechar='|')
    
    for sing_obs in singleton_table: #for each singleton observation
        if sing_obs == header_list: #don't waste time looping thru header a bunch
            keep_dict[keep_count] = header_list
            keep_count += 1
            continue #skip to next iteration
        else:
            for mask_chrom in masks_for_chrom.keys(): #loop through mask table 
                chrom = sing_obs[0]
                loc = int(sing_obs[1])
                sing_doub = sing_obs[2]
                the_allele = sing_obs[3]
                indv = sing_obs[4]
                keep_obs_list = [chrom, str(loc), sing_doub, the_allele, indv] #keep entire row/observation
                if chrom == 'LT635626' or  chrom == 'LT635627': #keep all apicoplast and mitochondrial observations
                    keep_dict[keep_count] = keep_obs_list
                    keep_count += 1
                elif chrom not in masks_for_chrom:  # don't crash if there are not masks for chrom
                    continue 
                else:
                    for start,end in masks_for_chrom[chrom]: #loop through mask table
                        if start < loc and loc < end:
                            if keep_obs_list not in keep_dict.values():
                                keep_dict[keep_count] = keep_obs_list
                                keep_count += 1
                            else:
                                pass
                        else:
                            pass
                

out_file = 'trimmed-' + singfile + '.csv'

with open(out_file, 'w') as out_file:
    out_writer = csv.writer(out_file, delimiter=",")

    for lineno, obs_list in keep_dict.items():
        out_writer.writerow(obs_list)
