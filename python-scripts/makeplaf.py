#!/usr/bin/env python3
#script modified from Vir's oroiginal script (/data/wraycompute/vdp5/scripts/deploid/makeplaf.py)

#pip install PyVCF to import vcf
import os
import vcf
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', help="Tab file for use")
args = parser.parse_args()


vcf_reader = vcf.Reader(open(args.vcf, 'r'))


newnombre = '.'.join(args.vcf.split('.')[:-1] + ['plaf',])
print(newnombre)

newfle = open(newnombre, 'w')
newfle.write('CHROM\tPOS\tMAF\n')
i = 0
for record in vcf_reader:
	if record.num_unknown == len(record.samples): #if all samples have missing genotypes, skip the position
		pass
	else:
#		newfle.write('{}\t{}\t{}\n'.format(record.CHROM, record.POS, record.aaf[0]))
                newfle.write(f'{record.CHROM}\t{record.POS}\t{record.aaf[0]}\n') #.format(record.CHROM, record.POS, record.aaf[0]))

newfle.close()

