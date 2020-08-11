from sys import argv
import re

with open(argv[1], 'r') as infile:
    print '\t'.join(["id","chr1", "pos1", "chr2", "pos2", "ori", "type","size","illumina","nanopore","pacbio","tenx", "pcr", "capture", "bionano"])
    for line in infile:
        if line.startswith('#'):
            continue
        t = "."
        l = "."
        ori = '.'
        ill = "0"
        ont = "0"
        pb = "0"
        tx = "0"
	pcr = "0"
	capture = "0"
	bn = "0"
        chr1, pos1, ident, ref, alt, qual, filt, info, fmt, gt = line.rstrip().split('\t')  
        for i in info.split(';'):
            if "SVLEN" in i:
                l = i.split('=')[-1]
            elif "SVTYPE" in i:
                t = i.split('=')[-1]
            elif "STRANDS" in i:
                ori = i.split('=')[-1]
        if alt != '<INS>':
            alt_match = re.search("^(\w*\]|\w*\[)(\w+):(\d+)(\]\w*|\[\w*)$", alt)
            chr2 = alt_match.group(2)
            pos2 = alt_match.group(3)
        else:
            chr2 = chr1
            pos2 = pos1
            
        if "ILLUMINA" in filt:
            ill = '1'
        if "NANOPORE" in filt:
            ont = '1'
        if "PACBIO" in filt:
            pb = '1'
        if "10X" in filt:
            tx = '1'
	if "PCR" in filt:
            pcr = '1'
	if "CAPTURE" in filt:
            capture = '1'
	if "BIONANO" in filt:
            bn = '1'

        print '\t'.join([ident,chr1, pos1, chr2, pos2, ori, t, l, ill, ont, pb, tx, pcr, capture, bn])
