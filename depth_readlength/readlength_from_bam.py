#! /bin/python

from sys import argv
import pysam


bamfile=argv[1]
chromosome=argv[2]
tech, sample=bamfile.split('.')[0:2]

samfile = pysam.AlignmentFile(bamfile, "rb")
idents = []
lengths = []

print "tech,sample,ident,length,chrom"
#Loops through all reads in chromosome and prints the length, skipping secondary alignments and duplicates through the read identifier
for read in samfile.fetch(chromosome):
    if read.is_secondary:
        continue
    ident = read.query_name
    if ident in idents:
        continue
    length = read.query_length
    lengths.append(length)
    idents.append(ident)
    print ','.join([tech, sample, ident, str(length), str(read.reference_name)])
