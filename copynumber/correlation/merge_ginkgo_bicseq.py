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


illumina_file = "../circos/illumina.tumor.circos"
nanopore_file = "../circos/nanopore.tumor.circos"
pacbio_file = "../circos/pacbio.tumor.circos"
tenx_file = "../circos/tenx.tumor.circos"
bionano_file = "../circos/bionano.tumor.circos"

masterD = {}
#use nanopore as template to start the dictionary
with open(nanopore_file, 'r') as infile:
    for line in infile:
        chrm, start, end, cn = line.rstrip().split("\t")
        if chrm == "hsX" or chrm == "hsY":
            continue
        masterD.setdefault(chrm, []).append([int(start), int(end), int(cn)])
for filename in [illumina_file,pacbio_file,tenx_file,bionano_file]:
    masterD = add_file_to_masterD(filename, masterD)

print '\t'.join(["chr","start","end","ONT","ILL","PB","TENX","BNG"])
for chrm, l in masterD.items():
    for i in l:
        print '\t'.join([chrm]+[str(x) for x in i])
