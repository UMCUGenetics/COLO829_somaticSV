cn_file = "COLO-829--COLO-829BL.cnv.annotated.v6.final.bed"

d = {}
with open(cn_file, 'r') as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        line = line.rstrip().split("\t")
        chrm, start,end = line[:3]
        ratio = line[4]
        #new ratio is converted to ploidy 3 and normalized to a diploid genome and to 0. So a copy number of 3 would be a 1.
        new_ratio = ((2**float(ratio)) * 3) - 2
        for pos in range(int(start), int(end)+1, 10000):
            print '\t'.join([chrm, str(pos), str(pos + 9999), str(new_ratio)])
