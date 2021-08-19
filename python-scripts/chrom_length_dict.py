
#!/usr/bin/env python3

""" Haploid Summary Statistics
Function to extract chromosome and length from your VCF and store information in a dictionary.
Returns dictionary in the format:
{'LT6356##': length}
"""

import re
import collections #default dict
import sys

def chrom_length_dict(vcf_file):
    chrom_len_dict = collections.defaultdict(str)
    with open(vcf_file, "r") as file:
        for line in file:
            if re.match(r"##contig=<ID=LT6356.+,", line): #keep only PvP01 chromosome lines
                #parse the line to get chrom and length
                chrom_name = re.search(r"LT6356[0-9][0-9]", line)
                chrom_length = re.search(r"length=[0-9]*", line) 
                
                #save to dictionary
                chrom_len_dict[chrom_name.group(0)] = chrom_length.group(0).split("=")[1] 
                """
                Notes:
                - print(chrom_name.group(0)) <- group(0) is how you access the matched part of the string 
                https://stackoverflow.com/a/41220180/10176950
                - for the value: split the length string by the '=' and only keep the second part of the split using [1]
                """
            else:
                pass # don't keep contigs or non-PvP01 chromosomes
    return chrom_len_dict


def pv_chrom_length_lookup(my_chrom):
    chrom_dict = {'LT635612': 1021664, 
    'LT635613': 956327, 
    'LT635614': 896704, 
    'LT635615': 1012024, 
    'LT635616': 1524814, 
    'LT635617': 1042791,
    'LT635618': 1652210,
    'LT635619': 1761288,
    'LT635620': 2237066,
    'LT635621': 1548844,
    'LT635622': 2131221,
    'LT635623': 3182763,
    'LT635624': 2093556,
    'LT635625': 3153402,
    'LT635626': 29582, # apicoplast
    'LT635627': 5989 # mitochondria
    }

    return chrom_dict[my_chrom]
  


if __name__=="__main__":
   input_file = sys.argv[1] 
   my_dict = chrom_length_dict(input_file)
   print(my_dict)