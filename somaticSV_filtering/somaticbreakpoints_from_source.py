#! /bin/sh
from sys import argv
import re


soms = argv[1]
headerfile = "somatic.merged.clean.vcf"
ill = "illumina.somatic.vcf"
pb = "pacbio.somatic.vcf"
ont = "nanopore.somatic.vcf"
tenx = "tenx.somatic.vcf"

illD = {}
tenxD = {}
pbD = {}
ontD = {}

def construct_alt(chr2, pos2, ori, ref):
    if ori == '++':
        alt = "%s]%s:%s]" % (ref, chr2, pos2)
    elif ori == '--':
        alt = "[%s:%s[%s" % (chr2, pos2, ref)    
    elif ori == '+-':
        alt = "%s[%s:%s[" % (ref, chr2, pos2) 
    elif ori == '-+':
        alt = "]%s:%s]%s" % (chr2, pos2, ref)
    return alt



with open(soms, 'r') as infile:
    for line in infile:
        line = line.rstrip().split('\t')
        if line[0] == 'svid':
            continue
        ident, chr1, pos1, chr2, pos2, svtype, result, ori, svlen, i, p, t, o, pcr, capt, bn, notrelevant, supportset = line
        if i == "1":
            illD[ident] = {'chr1':chr1, 'pos1':pos1, 'chr2':chr2, 'pos2':pos2, 'illumina':i, 'pacbio':p, 'tenx':t, 'nanopore':o, 'pcr':pcr, 'capture':capt, 'bionano':bn, 'svtype':svtype, 'supportset':supportset, 'result':result}
        elif t == "1":
            tenxD[ident] = {'chr1':chr1, 'pos1':pos1, 'chr2':chr2, 'pos2':pos2, 'illumina':i, 'pacbio':p, 'tenx':t, 'nanopore':o, 'pcr':pcr, 'capture':capt, 'bionano':bn, 'svtype':svtype, 'supportset':supportset, 'result':result}
        elif o == "1":
            ontD[ident] = {'chr1':chr1, 'pos1':pos1, 'chr2':chr2, 'pos2':pos2, 'illumina':i, 'pacbio':p, 'tenx':t, 'nanopore':o, 'pcr':pcr, 'capture':capt, 'bionano':bn, 'svtype':svtype, 'supportset':supportset, 'result':result}
        elif p == "1":
            pbD[ident] = {'chr1':chr1, 'pos1':pos1, 'chr2':chr2, 'pos2':pos2, 'illumina':i, 'pacbio':p, 'tenx':t, 'nanopore':o, 'pcr':pcr, 'capture':capt, 'bionano':bn, 'svtype':svtype, 'supportset':supportset, 'result':result}


with open(headerfile, 'r') as vcf:
    for line in vcf:
        if line.startswith("#CHROM"):
            print line.rstrip()
            break
        elif line.startswith("#"):
            print line.rstrip()
        else: 
            continue


def print_vcf(dictionary, vcffile):    
    with open(vcffile, 'r') as vcf:
        for l in vcf:
            if l.startswith('#'):
                continue
            l = l.rstrip().split('\t')
            chr1, pos1 = l[:2]
            ref = l[3]
            alt = l[4]
            info = l[7]
            # create regex match varianbles
            alt_match = re.search("^(\w*\]|\w*\[)(\w+):(\d+)(\]\w*|\[\w*)$", alt)
            chr2_match = re.search("CHR2=(.+?)(;|$)", info)
            end_match = re.search("END=(\d+)(;|$)", info)
            # set second chromosome
            
            chr2= chr1
            if chr2_match:
                    chr2 = chr2_match.group(1)
            elif alt_match:
                    chr2 = alt_match.group(2)
            chr2 = chr2.replace("chr","")
            
            # set end position
            if end_match:
                    pos2 = end_match.group(1)
            elif alt_match:
                    pos2 = alt_match.group(3)
                    
            if not alt_match:
                ori_match = re.search("STRANDS=([+-][+-](;|$))", info)
                if ori_match:   
                    ori = ori_match.group(1)
                    new_alt = construct_alt(chr2, pos2, ori, ref)
                else: #this is for PBSV vcf, which has the orientations all messed up
                    svtype_match=re.search("SVTYPE=(.+?)(;|$)", info)
                    svtype = svtype_match.group(1)
                    if svtype == "BND":
                        new_alt = alt
                    elif svtype == "INS":
                        new_alt = "<INS>"
                    elif svtype == "DEL":
                        ori = "+-"
                        new_alt = construct_alt(chr2, pos2, ori, ref)
                    elif svtype == "DUP":
                        ori = "-+"
                        new_alt = construct_alt(chr2, pos2, ori, ref)                    
                l[4] = new_alt
            
            for key, d in dictionary.items():
                if chr1 == d["chr1"]:
                    if chr2 == d["chr2"] or d['svtype'] == "INS":
                        if int(pos1) in range(int(d["pos1"])-50, int(d["pos1"])+51):
                            if int(pos2) in range(int(d["pos2"])-50, int(d["pos2"])+51) or d['svtype'] == "INS":
                                validations = []
                                if d['pcr'] == '1':
                                    validations.append('PCR')
                                if d['capture'] == '1':
                                    validations.append('CAPTURE')
                                if d['bionano'] == '1':
                                    validations.append('BIONANO')
                                if len(validations) == 0:
                                    validations=["NOTVALIDATED"]
                                l[7] += ";VALIDATIONS="+';'.join(validations)
                                l[7] += ";SUPPORTSET="+d['supportset']
                                l[6] = d['result']
                                print '\t'.join(l[:9] + ["0/1"])
                                del dictionary[key]
                                break
                    
print_vcf(illD, ill)
print_vcf(tenxD, tenx)
print_vcf(ontD, ont)
print_vcf(pbD, pb)



                        
    
