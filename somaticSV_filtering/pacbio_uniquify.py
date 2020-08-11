#! /bin/python

filename = "pacbio.somatic.adjusted.vcf"

with open(filename, 'r') as infile:
    mate_l = []
    for line in infile:
        line = line.rstrip()
        if line.startswith("#"):
            print line
            continue
        if not "MATEID" in line:
            print line
            continue
        ident = line.split('\t')[2]
        if ident in mate_l:
            continue
        mate = line.split('\t')[7].split('=')[-1]
        print line
        mate_l.append(mate)
