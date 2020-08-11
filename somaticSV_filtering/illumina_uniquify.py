#! /bin/python

filename = "illumina.somatic.adjusted.vcf"

with open(filename, 'r') as infile:
    l = []
    for line in infile:
        line = line.rstrip()
        if line.startswith("#"):
            print line
            continue
        ident = line.split('\t')[2][:-1]
        if ident in l:
            continue
        else:
            print line
            l.append(ident)
        