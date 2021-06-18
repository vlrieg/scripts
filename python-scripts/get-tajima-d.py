#!/usr/bin/env python3

import argparse
import math

# Tajima's D equation 5.5 from:
# https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
# running this script with values from the first example gives me the correct output:
# $ ./get-tajima-d.py --pi 3.888889 --segregating_sites 16 --no_samples 10
# -1.4461719856148525


# Important note on Hartl vs Tajima definitions of Pi and S (from link above):
# Hartl defines pi as the fraction of nt site differences and S as the fraction of sites segregating (SNPs/total bp)
# Tajima defines pi as the mean number of sites that differ and S as the number of sites segregating
# As long as both pi and S were calculated using EITHER Hartl's OR Tajima's method, the answer from this equation will be the same

parser = argparse.ArgumentParser(description="calculate Tajima's D")

parser.add_argument("--pi", help="value of pi") # args.pi
parser.add_argument("--segregating_sites", help="value for S") # args.segregating_sites
parser.add_argument("--no_samples", help="how many sequences/individuals you have") # args.no_samples

args = parser.parse_args()

### convert argparse input to variables
pi = float(args.pi)
s_sites = float(args.segregating_sites)
n_samps = int(args.no_samples)

### functions
def take_sum(n, exponent):
    my_sum = 0
    for sample in range(1,n):
        #print(sample) #make sure my range is from 1 to n-1
        my_sum += 1 / (sample ** exponent)
    return(my_sum)


### set up values
a1 = take_sum(n_samps, 1)
a2 = take_sum(n_samps, 2)

b1 = (n_samps + 1) / (3 * (n_samps - 1))
b2 = (2 * ((n_samps ** 2) + n_samps + 3)) / (9 * n_samps * (n_samps - 1))

c1 = b1 - (1 / a1)
c2 = b2 - ((n_samps + 2) / (a1 * n_samps)) + (a2 / (a1 ** 2))

e1 = c1 / a1
e2 = c2 / ((a1 ** 2) + a2)

### calculate Tajima's D and print
tajimad = (pi - (s_sites / a1)) / math.sqrt((e1 * s_sites) + ((e2 * s_sites) * (s_sites - 1)))

print(tajimad)