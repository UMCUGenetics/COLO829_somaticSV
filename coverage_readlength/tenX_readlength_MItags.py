#! /bin/python 

from sys import argv
import pysam

bamfile=argv[1]
chromosome=argv[2]
tech, sample=bamfile.split('.')[0:2]
samfile = pysam.AlignmentFile(bamfile, "rb")
print "tech,sample,ident,length,chrom"

mi_d = {}
#structure: mi_d = {tag1:{chr1:{start = [], end = []}, chr2: {start = [], end = []}, tag2 ...}}
for read in samfile.fetch(chromosome):
    if read.has_tag('MI'):
        mi = read.get_tag('MI')
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        if mi in mi_d:
            if chrom in mi_d[mi]:                        
                mi_d[mi][chrom]['starts'].append(start)
                mi_d[mi][chrom]['ends'].append(end)

            else:
                mi_d[mi][chrom] = {'starts':[start], 'ends':[end]}
        else:
            mi_d[mi] = {chrom:{'starts':[start], 'ends':[end]}}
                

for tag, d in mi_d.items():
    if len(d.keys()) > 1:
        ##Mapping to multiple chromosomes, discard this molecule
        continue
    else:
        chrom = d.keys()[0]
        beg = min(d[chrom]['starts'])
        end = max(d[chrom]['ends'])
        l = end - beg
        print ','.join([tech, sample, str(tag), str(l), str(chrom)])
        
