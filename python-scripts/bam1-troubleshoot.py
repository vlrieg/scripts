#!/usr/bin/env python
# Example command line: 
# ./bam1.py ~/mnt/hardac/gpfs/fs1/data/wraycompute/vdp5/reference_data/PVP01.fasta ~/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/ERR2679006_dedup_reads.bam LT635626

#refname = sys.argv[1] #reference sequence (i.e. PvP01)
#bamfilename = sys.argv[formerly2] #don't need this function argument after refactoring
#seqname = sys.argv[2] #which chromosome are you interested in?

import pysam
import sys
from collections import defaultdict
import operator
import csv
import glob #Unix style pathname pattern expansion
import re #regular expressions

#############################################################################################################################

#running scripts from local mount
#filenames = glob.glob('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax-like/bam/dedup/*_dedup_reads.bam')
#filenames = glob.glob('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/*_dedup_reads.bam')

#running scripts from hardac
filenames = glob.glob('/gpfs/fs1/data/wraycompute/malaria/vivax-like/bam/dedup/*_dedup_reads.bam')
#filenames = glob.glob('/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/*_dedup_reads.bam')

#pull out species name for output file
#try:
    #local mnt
    #species = re.search('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/(.+?)/(.+?)/(.+?).bam', filenames[0]).group(1)

    #on hardac
    #species = re.search('/gpfs/fs1/data/wraycompute/malaria/(.+?)/(.+?)/(.+?).bam', filenames[0]).group(1)
    #species = re.search('~/wraycompute/malaria/(.+?)/(.+?)/(.+?).bam', filenames[0]).group(1)
    #species = "vivax"
species = "vivax-like"

#except AttributeError:
#    print(f'species = {species}')

#############################################################################################################################

#Open CSV file to write output
file = open(f'{species}_{sys.argv[2]}_allele-freq_output.csv', mode='a')             
writer = csv.writer(file)
    
fieldnames = ['sample', 'loci', 'index', 'coverage', 'allele', 'frequency'] # column header names
writer.writerow(fieldnames)

#############################################################################################################################

def allele_freq(refname, bamfilename, seqname):
    bamfile = pysam.AlignmentFile(bamfilename,"rb")

    # #pull out sample name for output file
    try:    
	# #local hardac mnt
	# #SRArun = re.search('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax-like/bam/dedup/(.+?).bam', filenames[0]).group(1)
    # #SRArun = re.search('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/(.+?).bam', filenames[0]).group(1)
    # #on hardac
    	SRArun = re.search('/gpfs/fs1/data/wraycompute/malaria/vivax-like/bam/dedup/(.*).bam', bamfilename)
        #search_obj = re.search('/mnt/hardac/gpfs/fs1/data/wraycompute/malaria/vivax-like/bam/(.*)/(.*).bam', bamfilename)
        #SRArun = search_obj.group(2)
    #     SRArun = re.search('/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/(.+?).bam', filenames[0]).group(1)
    except AttributeError:
        print(f'SRArun = {SRArun}')

    #SRArun = 'ERR2679006_dedup_reads'

    #calculate major allele info
    sites = 0
    for pileupcolumn in bamfile.pileup(contig=seqname):
        sites = sites + 1
        counts = defaultdict(int)
        total = 0
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                letter = pileupread.alignment.query_sequence[pileupread.query_position]
                counts[letter] = counts[letter] + 1
                total = total + 1
        if len(counts) < 2:
            continue

        major_allele = max(counts.items(), key=operator.itemgetter(1))[0]
        major_allele_count = counts[major_allele]
        major_allele_fraction = major_allele_count/total

        if major_allele_fraction < 0.99 and major_allele_count < total - 10:
            results = [SRArun, seqname, pileupcolumn.pos, pileupcolumn.n, major_allele, major_allele_fraction]
            writer.writerow(results)

    print(f'{bamfilename}')
    #print(f"scanned {sites} sites for {SRArun}")
    bamfile.close()

#############################################################################################################################

for x in filenames:
    allele_freq(sys.argv[1], x, sys.argv[2])

#file_006 = '/gpfs/fs1/data/wraycompute/malaria/vivax/DuffyNeg/ERR2679006_dedup_reads.bam'
#allele_freq(sys.argv[1], file_006, sys.argv[2])

file.close()
