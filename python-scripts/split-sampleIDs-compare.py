#!/usr/bin/env python3

# run like:
# ./split-sampleIDs-compare.py ~/Dropbox\ \(Duke\ Bio_Ea\)/coalescent/samples/tables/sample_ID_accessions.csv > out.txt
# the sample IDs printed in 'out.txt' are those that are present in the data table twice 
# (one with run accession number ERR and another under the sample accession number ERS)
# delete the ERS rows that are duplicated in merged-sample-table.txt

# input file is sample_ID_accessions.csv generated by ~/scripts/r-scripts/find_duplicates_table_cleanup.Rmd
# input file is in the form 'accession,sample_ID'


import csv
import sys
from collections import defaultdict

acc_file = sys.argv[1]

with open(acc_file, "r") as infile:
    data = infile.read().splitlines() #this is a list

data_dict = defaultdict(str)

#accesion,sample_ID
for i in data:
    x = i.split(",")
    # print(x[0])
    # print(x[1])
    data_dict[x[0]] = x[1]

# for key, val in data_dict.items():
#    print(f"accession is {key} and sample ID is {val}")


ERR_list = []
ERS_list = []
NA_list = []

for accession in data_dict.keys():
    #print(accession[:3])
    if accession[:3] == 'ERR':
        ERR_list.append(accession)
    elif accession[:3] == 'ERS':
        ERS_list.append(accession)
    else:
        NA_list.append([accession, data_dict[accession]]) # list of lists where each list item is [accession, sample_ID]


# new dictionaries
err_dict = defaultdict(str)
for i in ERR_list:
    #print(i)
    #print(data_dict[i])
    err_dict[i] = data_dict[i]

ers_dict = defaultdict(str)
for i in ERS_list:
    ers_dict[i] = data_dict[i]


# check totals
#----
# err_dict_count = 0
# for key, val in err_dict.items():
#     err_dict_count += 1

# ers_dict_count = 0
# for key, val in ers_dict.items():
#     ers_dict_count += 1
# print(f'{len(ERR_list)} accessions in ERR list and {err_dict_count} accessions in the ERR dictionary')
# print(f'{len(ERS_list)} accessions in ERR list and {ers_dict_count} accessions in the ERS dictionary')
#----

this_count = 0
this_list = []
err_no_match_list = []
for i in ERR_list:
    #print(i)
    #print(err_dict[i])
    if err_dict[i] in ers_dict.values():
        this_count += 1
        #print(f'match: {err_dict[i]}')
        this_list.append(err_dict[i])
    else:
        err_no_match_list.append(err_dict[i])

# expected_remaineder = 162 - 96
# print(expected_remaineder) #66
# print(len(err_no_match_list)) #66


that_count = 0
that_list = []
for i in ERS_list:
    if ers_dict[i] in err_dict.values():
        that_count += 1
        #print(f'another match: {ers_dict[i]}')
        that_list.append(ers_dict[i])


this_list.sort()
that_list.sort()

for key, val in ers_dict.items():
    if val in that_list:
        print(key)

#############
# for i in range(0,96):
#     print(f'{this_list[i]},{that_list[i]}')




# more counting to double check

# match = 0
# no_match = 0
# for i in range(0,96):
#     if this_list[i] == that_list[i]:
#         match += 1
#         print('exact match')
#     else:
#         no_match += 1
#         print('no match')
    
#print(f'{match} matches and {no_match} mismatches')
#96 matches and no mismatches

#print(this_count)
#print(that_count)
#both counts are 96