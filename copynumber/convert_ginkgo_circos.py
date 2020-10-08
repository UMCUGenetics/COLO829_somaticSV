#! /bin/python

import sys

input = sys.argv[1]
column = int(sys.argv[2])
with( open( input, 'r' ) ) as f:
	f.next()
        for line in f:
            line=line.rstrip().split('\t')
            chrom, start, end = line[:3]
            #Normalize to 0
            cn = int(line[column]) - 2
            if cn < 0:
                cn = 0
            print("hs"+chrom.replace("chr","")+"\t"+str(start)+"\t"+str(end)+"\t"+str(cn))
