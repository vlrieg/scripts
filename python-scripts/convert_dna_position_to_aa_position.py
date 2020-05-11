#!/usr/bin/env python3

import sys
import math

#convert DNA position to AA position by dividing DNA position by 3 and then round up to next integer

dna_position = sys.argv[1]

aa_position = math.ceil(float(dna_position)/3)
print(f'DNA position {dna_position}\nis equeal to\nAmino Acid postiion {aa_position}')