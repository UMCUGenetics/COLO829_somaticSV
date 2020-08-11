from __future__ import division



def add_file_to_masterD(filename,d):
    newD = {}
    with open(filename, 'r') as infile:
        for line in infile:
            chrm, start, end, cn = line.rstrip().split("\t")
            newD.setdefault(chrm, []).append([int(start), int(end), int(round(float(cn)))])

    for chrm, l in d.items():
        for i, segment in enumerate(l):
            start, end = segment[:2]
            candidates = []
            for segment2 in newD[chrm]:
                start2, end2, cn2 = segment2
                if start2 > end or end2 < start:
                    continue
                candidates.append(cn2)
            if len(candidates) == 0:
                new_cn = ""
            else:
                new_cn = int(round(sum(candidates)/len(candidates)))
            d[chrm][i].append(new_cn)

    #Fill the empty ones with the consecutive cn
    for chrm, l in d.items():
        for i, segment in enumerate(l):
            if segment[-1] == "":
                #Pick the next one
                for segment2 in l[i+1:]:
                    if segment2[-1] != "":
                        segment[-1] = segment2[-1]
                        break
                #If still empty, pick the previous
                for segment2 in l[i-1::-1]:
                    if segment2[-1] != "":
                        segment[-1] = segment2[-1]
                        break

    return d


illumina_file = "illumina.tumor.circos"
nygc_file = "NYGC_COLO829.circos"
scA_file = "singleCell_groupA_cnv.circos"
scB_file = "singleCell_groupB_cnv.circos"
scC_file = "singleCell_groupC_cnv.circos"
scD_file = "singleCell_groupD_cnv.circos"
nanopore_file="nanopore.tumor.circos" #To be used as guide for the bins

masterD = {}
#use nanopore as template to start the dictionary
with open(nanopore_file, 'r') as infile: #Only as BIN guide, not interested in this copy number
    for line in infile:
        chrm, start, end, cn = line.rstrip().split("\t")
        if chrm == "hsX" or chrm == "hsY":
            continue
        masterD.setdefault(chrm, []).append([int(start), int(end)])
for filename in [illumina_file, nygc_file, scA_file, scB_file, scC_file, scD_file]:
    masterD = add_file_to_masterD(filename, masterD)

print '\t'.join(["chr","start","end","ILL","NYGC","scA","scB","scC","scD"])
for chrm, l in masterD.items():
    for i in l:
        print '\t'.join([chrm]+[str(x) for x in i])
