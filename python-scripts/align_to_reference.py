#!/usr/bin/env python3

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

alignment_filename = sys.argv[1]
contigs_filename = sys.argv[2]


# 1. Load all the contigs into a dictionary
seq = {}
with open(contigs_filename, 'rb') as fastafile:
    for record in SeqIO.parse(fastafile, "fasta"):
        seq[record.id] = record

#2. Look all alignments
orientations = {}
fractions = {}
start_point = {}
end_point = {}
alignment = {}
for multiple_alignment in AlignIO.parse(alignment_filename, "maf"):
#    print("\nprinting a new multiple alignment")

    chr_part = multiple_alignment[0]
    cont_part = multiple_alignment[1]

    chr_name = chr_part.id
    chr_match_length = int(chr_part.annotations["size"])
    chr_start = int(chr_part.annotations["start"])
    chr_end = chr_start + chr_match_length - 1
    chr_length = int(chr_part.annotations["srcSize"])

    contig_name = cont_part.id
    contig_match_length = int(cont_part.annotations["size"])
    contig_start = int(cont_part.annotations["start"])
    orientation = cont_part.annotations["strand"]
    contig_end = contig_start + orientation * (contig_match_length-1)
    contig_length = int(cont_part.annotations["srcSize"])
#    print("match_length = {}  length = {}".format(contig_match_length,contig_length))

    alignment[chr_name] = multiple_alignment

    if chr_name not in start_point:
        orientations[chr_name] = {}
        fractions[chr_name] = {}
        start_point[chr_name] = {}
        end_point[chr_name] = {}

    if contig_name not in start_point[chr_name]:
        orientations[chr_name][contig_name] = orientation
        start_point[chr_name][contig_name] = chr_start
        end_point[chr_name][contig_name] = chr_end
        fractions[chr_name][contig_name] = 0

    fractions[chr_name][contig_name] += float(contig_match_length)/float(contig_length)
    start_point[chr_name][contig_name] = min(start_point[chr_name][contig_name], chr_start)
    end_point[chr_name][contig_name] = max(end_point[chr_name][contig_name], chr_end)
#    print("fraction {} of contig {} now maps to {}  (orientation = {})".format(fractions[chr_name][contig_name],contig_name, chr_name,orientation))


# 3. Concatenate the sequences
for chr_name in start_point:
    lastpoint = None
    sequence = ""
    description = ""
    assembled_length = 0
    for contig_name,start in sorted(start_point[chr_name].iteritems(),key=lambda (k,v): (v,k)):
        fraction = fractions[chr_name][contig_name]
        useq = seq[contig_name].seq
        orientation = "+" if (orientations[chr_name][contig_name]==1) else "-";
        fraction = fractions[chr_name][contig_name]
        start = start_point[chr_name][contig_name]
        end = end_point[chr_name][contig_name]
#        print "{}/{}:{}   {}".format(chr_name, contig_name, start, fraction)

        if fraction < 0.9:
            continue
        if orientation == "-":
            useq = useq.reverse_complement()
        elif orientation == "+":
            pass
        else:
            raise

        if lastpoint is not None:
            if end <= lastpoint:
                raise "This shouldn't happen!"
            if start <= lastpoint:
                overlap = lastpoint - start + 1
                useq = useq[overlap:]
                start += overlap
 #               print "{}/{}: overlap = {}".format(chr_name,contig_name,overlap)
            
        sequence += str(useq)
        description += "  {} [{} - {} ({})]".format(contig_name, start, end, orientation)
        assembled_length += (end - start + 1)
        lastpoint = end

    reference_length = alignment[chr_name][0].annotations["srcSize"]

    description = "[assembled {}/{}={:.2f}%  missing {}] {}".format(assembled_length,
                                                            reference_length,
                                                            100*float(assembled_length)/
                                                                    float(reference_length),
                                                               reference_length-assembled_length,
                                                               description)
    
    record = SeqRecord(Seq(sequence), id=chr_name, name=chr_name, description=description)
    with open(chr_name + ".fasta",'w+') as outfile:
        SeqIO.write(record, outfile, "fasta")
