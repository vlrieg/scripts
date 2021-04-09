#!/usr/bin/env python3

# This script is meant to be run with a file of aligned sequences in fasta format
# e.g. P.vivax.mtDNA.muscle.output
# first argument is the input file

from Bio import SeqIO #to parse the fasta file
import sys # used so the file name can be entered on the command line as an argument
import datetime #to keep track of run times
import math # used to calculated standard deviation

# histogram packages
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


# record date and time when script started running
d_date = datetime.datetime.now()
start_time = d_date.strftime("%Y-%m-%d %I:%M:%S %p")

# SeqIO is a generator that imports fasta records as a list of seq objects with 
# properties "id" and "sequence"
mtDNA = list(SeqIO.parse(sys.argv[1], "fasta")) # sys.argv[1] is the name of the fasta file, w/o quotes
mtDNA_length = len(mtDNA)

pi_list = [] # counter to add up all the pi values
loop_count = 0 # counter to keep track of the number of sequence comparisons

for index_1 in range(mtDNA_length):           # for loop to index the first sequence
    for index_2 in range(mtDNA_length):   # for loop to index the second sequence
        if index_1 > index_2:                 # excludes index_1 = index_2 AND any duplicate comparisons
            #print(index_1, index_2)
            seq_1 = mtDNA[index_1]            # rename the sequences for simplicity
            seq_2 = mtDNA[index_2]            

            # ensure that sequence 1 is the same length as sequence 2
            # the argument is the message displayed if the assertion fails
            assert len(seq_1) == len(seq_2), f"""The sequences are NOT the same length!
            {seq_1} is {len(seq_1)} bases long and {seq_2} is {len(seq_2)} bases long"""
            
            # counters
            unknown_count = 0
            agree_count = 0
            disagree_count = 0
            
            # count the # of agrees, disagrees, and unknowns per base
            # for the two sequences
            for i in range(len(seq_1)):
                if seq_1[i] == "-" or seq_2[i] == "-":
                    unknown_count += 1
                elif seq_1[i] == "N" or seq_2[i] == "N":
                    unknown_count += 1
                elif seq_1[i] == seq_2[i]:
                    agree_count += 1
                else:
                    disagree_count += 1
            total_count = unknown_count + agree_count + disagree_count

            # ensure every base was counted in the loop
            assert total_count == len(seq_1), f"""Something went wrong! 
            total_count = {total_count} and sequence 1 length = {len(seq_1)} for 
            {mtDNA[index_1]} and {mtDNA[index_2]}"""

            # calculate pi for the current sequence comparison
            pi = disagree_count / (agree_count + disagree_count)

            # add up all values of pi (for all seqs)
            pi_list.append(pi)
            loop_count += 1
    print(f"finished sequence {index_1}")


# calculate pi
# disagree / agree + disagree
# then average over the number of comparisons
# 0.001-0.003 is reasonable (0.001 is more reasonable)
pi_sum = sum(pi_list)
pi_avg = pi_sum / loop_count  
print(f"The average value of pi for {loop_count} comparisons of {mtDNA_length} sequences is:", pi_avg)

# calculate the standard deviation
variance_list = [] #empty list for each instance of: pi - the ave pi value, squared

for pi in range(len(pi_list)):
    mean_dist_sq = (pi_list[pi] - pi_avg) ** 2
    variance_list.append(mean_dist_sq)

var_sum = sum(variance_list)
standard_deviation = math.sqrt(var_sum / loop_count)

print(f"The standard deviation is: {standard_deviation}")

#record date and time when script finished running
e_date = datetime.datetime.now()
end_time = e_date.strftime("%Y-%m-%d %I:%M:%S %p")

# script duration
print("Script started:", start_time)
print("Script ended: ", end_time)

# make a histogram for values of pi
num_bins = 20
n, bins, patches = plt.hist(pi_list, num_bins, facecolor = 'blue', alpha = 0.5)

plt.xlabel('Pi Value')
plt.ylabel('Count')
plt.title("Histogram of Pi values for global sequence comparisons")

plt.savefig("global-pi-value-histogram.pdf")
