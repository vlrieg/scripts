#!/usr/bin/env python3

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]


# Read the two files into lists
with open(file1, 'rb') as original_file:
    lines1 = original_file.read().decode('utf-8').splitlines()

with open(file2, 'rb') as new_file:
    lines2 = new_file.read().decode('utf-8').splitlines()


# are there duplicates within each list?
nodupes1 = list(dict.fromkeys(lines1))
nodupes2 = list(dict.fromkeys(lines2))

original1_len = len(lines1)
nodupe1_len = len(nodupes1)

original2_len = len(lines2)
nodupe2_len = len(nodupes2)

print(f'file 1 ({file1}) original length: {original1_len} accession numbers vs no duplicates length: {nodupe1_len}')
print(f'file 2 ({file2}) original length: {original2_len} acession numbers vs no duplicates length: {nodupe2_len}')

import collections
# print(f'duplicate samples in {file1}:')
# print([item for item, count in collections.Counter(lines1).items() if count > 1])
# print(f'duplicate samples in {file2}:')
# print([item for item, count in collections.Counter(lines2).items() if count > 1])


def comparison(list1, list2, unique_or_dup):
        accession_list = []
        if unique_or_dup == 'unique': #finds accessions that are unique to list1
                for accession in list1:
                        if accession not in list2:
                                accession_list.append(accession)
        elif unique_or_dup == 'duplicate': #finds accessions in list1 that are also in list2
                for accession in list1:
                        if accession in list2:
                                accession_list.append(accession)
        return accession_list

#print(f'Accessions unique to {file1} are:')
# for accession in unique_accession1:
#     print(accession)
# print()
unique_1 = comparison(lines1, lines2, 'unique')
# for x in unique_1:
#         print(x)

print()
print()
# print(f'Accessions unique to {file2} are:')
# unique_2 = comparison(lines2, lines1, 'unique')
# for x in unique_2:
#         print(x)

# # Duplicates
# print(f'Accessions in {file1} that are ALSO in {file2} are:')
# duplicates1 = comparison(lines1, lines2, 'duplicate')
# for x in duplicates1:
#         print(x)
# print(f'Accessions in {file2} that are ALSO in {file1} are:')
# duplicates2 = comparison(lines2, lines1, 'duplicate')
# for x in duplicates2:
#         print(x)






print()
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print()

print()
print('UNIQUE') 

# #which are unique vs duplicated between lists?
unique_accession1 = []
unique_accession2 = []

for line in lines1:
    #print(line)
    if line not in lines2:
        #print(line)
        unique_accession1.append(line)

for line in lines2:
    if line not in lines1:
        unique_accession2.append(line)

length1 = len(unique_accession1)
length2 = len(unique_accession2)

print(f'file {file1} has {length1} accession numbers')
print(f'file {file2} has {length2} accession numbers')

print(f'{len(unique_accession1)} accessions unique to {file1}:')
for accession in unique_accession1:
    print(accession)

print()

print(f'{len(unique_accession2)} accessions unique to {file2}:')
for accession in unique_accession2:
    print(accession)

print()
print('DUPLICATES')    

duplicate_accessions1 = []
duplicate_accessions2 = []

for line in lines1:
    if line in lines2:
        duplicate_accessions1.append(line)

for line in lines2:
    if line in lines1:
        duplicate_accessions2.append(line)



print(f'{len(duplicate_accessions1)} accessions in {file1} are ALSO in {file2}:')
for accession in duplicate_accessions1:
   print(accession)
print()
print()


print(f'{len(duplicate_accessions2)} accessions in {file2} are ALSO in {file1}:')
for accession in duplicate_accessions2:
   print(accession)

    