#!/usr/bin/env python3

import sys
import math

#convert DNA position to AA position by dividing DNA position by 3 and then round up to next integer

dna_position = sys.argv[1]

def find_aa(pos):
    aa_position = math.ceil(float(pos)/3)
    return aa_position


output = find_aa(dna_position)

print(f'DNA position {dna_position}\nis equeal to\nAmino Acid postiion {output}')