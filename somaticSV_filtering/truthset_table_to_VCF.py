#! /bin/python

from sys import argv

def print_header():
    header = []
    header.append('##fileformat=VCFv4.1')
    header.append('##contig=<ID=1,length=249250621>')
    header.append('##contig=<ID=2,length=243199373>')
    header.append('##contig=<ID=3,length=198022430>')
    header.append('##contig=<ID=4,length=191154276>')
    header.append('##contig=<ID=5,length=180915260>')
    header.append('##contig=<ID=6,length=171115067>')
    header.append('##contig=<ID=7,length=159138663>')
    header.append('##contig=<ID=8,length=146364022>')
    header.append('##contig=<ID=9,length=141213431>')
    header.append('##contig=<ID=10,length=135534747>')
    header.append('##contig=<ID=11,length=135006516>')
    header.append('##contig=<ID=12,length=133851895>')
    header.append('##contig=<ID=13,length=115169878>')
    header.append('##contig=<ID=14,length=107349540>')
    header.append('##contig=<ID=15,length=102531392>')
    header.append('##contig=<ID=16,length=90354753>')
    header.append('##contig=<ID=17,length=81195210>')
    header.append('##contig=<ID=18,length=78077248>')
    header.append('##contig=<ID=19,length=59128983>')
    header.append('##contig=<ID=20,length=63025520>')
    header.append('##contig=<ID=21,length=48129895>')
    header.append('##contig=<ID=22,length=51304566>')
    header.append('##contig=<ID=X,length=155270560>')
    header.append('##contig=<ID=Y,length=59373566>')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append('##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV>"')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">')
    header.append('##INFO=<ID=SUPP_SEQ,Number=.,Type=String,Description="Vector of supporting technologies">')
    header.append('##INFO=<ID=SUPP_VAL,Number=.,Type=String,Description="Vector of supporting technologies">')
    header.append('##INFO=<ID=GENE,Number=.,Type=String,Description="Gene(s) overlap in breakpoint">')
    header.append('##INFO=<ID=CLUSTER,Number=.,Type=String,Description="IDs of other clustered breakpoints">')
    header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOLO829T')
    return header


def build_alt(chrm, pos, ref, ori):
    if ori == "++":
        alt = "%s]%s:%s]" % (ref, chrm, pos)
    elif ori == "--":
        alt = "[%s:%s[%s" % (chrm, pos, ref)
    elif ori == "+-":
        alt = "%s[%s:%s[" % (ref, chrm, pos)
    elif ori == "-+":
        alt = "]%s:%s]%s" % (chrm, pos, ref)
    return alt

def build_SUPP_SEQ(ill, ont, pb ,tenx):
    supp_seq = ""
    if ill == "1":
        supp_seq += "ILL,"
    if ont == "1":
        supp_seq += "ONT,"
    if pb == "1":
        supp_seq += "PB,"
    if tenx == "1":
        supp_seq += "10X,"
    supp_seq = supp_seq[:-1]
    return supp_seq


def build_SUPP_VAL(pcr, capture, bn):
    supp_val = ""
    if pcr == "1":
        supp_val += "PCR,"
    if capture == "1":
        supp_val += "CAPTURE,"
    if bn == "1":
        supp_val += "BIONANO,"
    if supp_val == "":
        supp_val = "NOT_VALIDATED,"
    supp_val = supp_val[:-1]
    return supp_val

def main():
    truthset_file = argv[1]
    line_list = []
    
    with open(truthset_file, 'r') as infile:
        infile.next()
        for line in infile:
            line=line.rstrip().split('\t')
            new_id, chr1, pos1, chr2, pos2, ref, ori, svtype, svsize, ill, ont, pb, tenx, pcr, capture, bn, gt, gene, cluster, comment = line
            if chr2 == "<INS>":
                alt1 = "<INS>"
            else:
                alt1 = build_alt(chr2, pos2, ref, ori)
                ori2 = ori[::-1]
                alt2 = build_alt(chr1, pos1, ref, ori2)
            new_id_1=new_id+"_1"
            new_id_2=new_id+"_2"
            supp_seq = build_SUPP_SEQ(ill, ont, pb, tenx)
            supp_val = build_SUPP_VAL(pcr, capture, bn)
            info = "SVLEN=%s;SVTYPE=%s;SUPP_SEQ=%s;SUPP_VAL=%s;GENE=%s;CLUSTER=%s" % (svsize, svtype, supp_seq, supp_val, gene, cluster.replace(" ",""))            
            record1 = [chr1, pos1, new_id_1, ref, alt1, ".", "PASS", info, "GT", gt]
            line_list.append(record1)
            if chr2 != "<INS>":
                record2 = [chr2, pos2, new_id_2, ref, alt2, ".", "PASS", info, "GT", gt]
                line_list.append(record2)
    header = print_header()
    for line in header: 
        print line
    for line in line_list:
        print  '\t'.join(line)
    



if __name__ == '__main__':
    main()
