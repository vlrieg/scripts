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

# Looks like:
# LT635612 [[0, 116540], [903590, 1021663]]
# LT635613 [[0, 100045], [745643, 956326]]
# LT635614 [[0, 108061], [894722, 896703]]
# LT635615 [[0, 185114], [967650, 1012023], [685685, 748923]]
# LT635616 [[0, 143101], [1408236, 1524813]]
# LT635617 [[0, 39348], [1017209, 1042790]]
# LT635618 [[0, 52470], [1463494, 1652209]]
# LT635619 [[0, 28423], [1627672, 1761287]]
# LT635620 [[0, 255761]]
# LT635622 [[0, 62171], [2051114, 2131220]]
# LT635623 [[0, 54908], [792292, 818496], [3028218, 3182762]]
# LT635624 [[0, 31856], [2038385, 2093555]]
# LT635625 [[0, 57], [3019711, 3153401]]



def in_masked_region(x, range_start, range_end):
    if int(x) < int(range_start) or int(x) > int(range_end):
        return False # x is NOT within masked region 
    else:
        return True # x is within masked region 


def keep_unmasked(obs_chrom, obs_pos):
    in_any_masked_region_list = []

    for chrom_key, masked_val in masks_for_chrom.items():
        if obs_chrom == chrom_key:
            for each_range in masked_val: #most chromosomes have more than one masked region
                if in_masked_region(obs_pos, each_range[0], each_range[1]) == False: #position NOT in this masked region
                    pass
                elif in_masked_region(obs_pos, each_range[0], each_range[1]) == True: #position is in this masked region
                    in_any_masked_region_list.append(1)
                else:
                    print(f"something weird is happening for observation {obs_chrom}: {obs_pos} queried against mask: {chrom_key}: {each_range}")
        else:
            pass #ignore the rest of the masks_for_chrom dictionary entries

    #want to check to see if there are more than one masked regions on a chrom before determining whether to keep or not
    if sum(in_any_masked_region_list) == 0: #the position is not in any masked region
        return True #keep this position
    else:
        return False # this position is in one of the masked regions for this chromosome


keep_dict = defaultdict(str)
keep_count = 0
header_list = ['CHROM', 'POS', 'SINGLETON/DOUBLETON', 'ALLELE', 'INDV']


# # Read singleton file into a list and then compare to mask table
# # only keep lines that don't fall within masked regions
# # also keep all lines that are from apicoplast and mitochondria

with open(singfile, mode='r') as singleton_file:
    singleton_table = csv.reader(singleton_file, delimiter = '\t', quotechar='|')
    
    for sing_obs in singleton_table: #for each singleton observation
        if sing_obs == header_list: #don't waste time looping thru header a bunch
            keep_dict[keep_count] = header_list
            keep_count += 1
            continue #skip to next iteration
        else:
            chrom = sing_obs[0]
            loc = int(sing_obs[1])
            sing_doub = sing_obs[2]
            the_allele = sing_obs[3]
            indv = sing_obs[4]
            keep_obs_list = [chrom, str(loc), sing_doub, the_allele, indv] #keep entire row/observation
            if chrom == 'LT635626' or  chrom == 'LT635627' or chrom == 'LT635621': #keep all apicoplast, mitochondrial, and chrom LT635621 observations
                if keep_obs_list not in keep_dict.values(): #only record the observation once
                    keep_dict[keep_count] = keep_obs_list
                    keep_count += 1
                else:
                    pass
            # elif chrom not in masks_for_chrom:  # don't crash if there are not masks for chrom
            #     continue 
                # not sure I need this elif statement if I'm accounting for all unmasked chroms above 
                # but I'm gonna leave it in just in case
            elif keep_unmasked(chrom, loc) == True:
                if keep_obs_list not in keep_dict.values(): #only record the observation once
                    keep_dict[keep_count] = keep_obs_list
                    keep_count += 1
                else:
                    pass
            else:
                pass
                # print()
                # print(f"in masked reange?:")
                # print(keep_obs_list)
                # print()


###to print values instead of saving to file:
# for key,val in keep_dict.items():
#     print(*val, sep="\t")

###to save to file:
out_file = 'trimmed-' + singfile + '.csv'

with open(out_file, 'w') as out_file:
    out_writer = csv.writer(out_file, delimiter=",")

    for lineno, obs_list in keep_dict.items():
        out_writer.writerow(obs_list)
