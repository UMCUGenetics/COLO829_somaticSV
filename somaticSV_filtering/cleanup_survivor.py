import vcf as pyvcf
import sys

def get_alt(record):
    chr2 = str(record.INFO['CHR2'])
    pos2 = str(record.INFO['END'])
    ref = str(record.REF[0])
    ori = str(record.INFO['STRANDS'])

    return chr2, pos2, ori
   
def get_support(support_vector):
    supp_seq_raw = support_vector[:5]
    supp_seq_raw = [ int(x) for x in supp_seq_raw]
    ill, ont, pb, tenx = supp_seq_raw
    supp_seq_list = []
    if ill:
        supp_seq_list.append('illumina')
    if pb:
        supp_seq_list.append('pacbio')
    if tenx:
        supp_seq_list.append('10X')
    if ont:
        supp_seq_list.append('nanopore')
    supp_seq = ','.join(supp_seq_list)
    return supp_seq
   
   
def inferr_svlen(record):
    svlen = 0
    if 'SVLEN' in record.INFO:
	if isinstance(record.INFO['SVLEN'], list):
            svlen = abs(int(record.INFO['SVLEN'][0]))
        else:
    	    svlen = abs(int(record.INFO['SVLEN']))
    if 'CHR2' in record.INFO:
        if record.CHROM == record.INFO['CHR2'] and svlen == 0:
            svlen = abs(int(record.INFO['END'])-int(record.POS))
    return svlen

    
def inferr_svtype(strands, type_raw, chr1, chr2):
    if type_raw == 'BND':
        t = 'BND'
    elif type_raw == 'INS':
        t = 'INS'
    else:
        if chr1 != chr2:
            t = 'BND'
        elif strands == '++':
            t = 'INV3'
        elif strands == '--':
            t = 'INV5'
        elif strands == '+-':
            t = 'DEL'
        elif strands == '-+':
            t = 'DUP'
    return t
   
def modify_info(record, ins_seq):
    supp=record.INFO['SUPP']
    supp_seq = get_support(record.INFO['SUPP_VEC'])
    end=record.INFO['END']
    chr2=record.INFO['CHR2']
    svlen = inferr_svlen(record)
    if ins_seq != '' and ins_seq != '<INS>':
        svlen = len(ins_seq)
    svtype_inferred = inferr_svtype(record.INFO['STRANDS'], record.INFO['SVTYPE'], record.CHROM, record.INFO['CHR2'])
    cipos = ','.join(record.INFO['CIPOS'])
    ciend = ','.join(record.INFO['CIEND'])
    strands = record.INFO['STRANDS']
    info = "SUPP=%s;SUPP_SEQ=%s;CHR2=%s;END=%s;SVLEN=%i;SVTYPE_INF=%s;CIPOS=%s;CIEND=%s;STRANDS=%s" % (supp,supp_seq, chr2, end, svlen, svtype_inferred, cipos, ciend, strands)
    if ins_seq != '':
        info += ";INS_SEQ=%s" % (ins_seq)
    return info

def extract_genotype(record):
    gt = record.samples[0]['GT']
    if gt != '1/1':
        gt = '0/1'
    return gt

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


def main():
    with open(sys.argv[1], 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                line= line.rstrip()
                flag = True
                if 'contig' in line and '_' in line:
                    flag = False
                elif 'contig' in line and 'ID=M,' in line:
                    flag = False
                elif 'contig' in line  and 'ID=hs37d5,' in line:
                    flag = False
                elif 'ALT' in line:
                    flag=False
                elif 'FORMAT' in line and not 'GT' in line:
                    flag=False
                ##CHECK INFO FIELDS  
                elif 'INFO' in line and '=MAPQ,' in line:
                    flag=False
                elif 'INFO' in line and '=RE,' in line:
                    flag=False
                elif 'INFO' in line and '=AVGLEN,' in line:
                    flag=False
                elif 'INFO' in line and '=SVMETHOD,' in line:
                    flag=False
                elif 'INFO' in line and '=SVTYPE,' in line:
                    flag=False
                elif 'INFO' in line and '=SUPP_VEC,' in line:
                    flag=False
                elif 'INFO' in line and '=SUPP,' in line:
                    flag=False
                
                elif line.startswith('#CHROM'):
                    flag = False
                if flag:
                    print line
                    
    print '##INFO=<ID=INS_SEQ,Number=.,Type=String,Description="Inserted sequence">'
    print '##INFO=<ID=SUPP_SEQ,Number=.,Type=String,Description="Vector of supporting technologies">'
    print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">'
    print '##INFO=<ID=SVTYPE_INF,Number=1,Type=String,Description="Type of the SV (simply inferred from orientation)">'
    print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOLO829T'
    
    vcf_reader = pyvcf.Reader(open(sys.argv[1], 'r'))
    new_records_d = {}
    counter = 1
    for record in vcf_reader:
        chr1 = str(record.CHROM)
        pos1 = str(record.POS)
        ident = str(record.ID)
        ref = str(record.REF[0])
        svtype = record.INFO['SVTYPE']
        ins_seq = ''
        alt = ''
        if svtype == "INS":
            alt = '<INS>'
            chr2=chr1
            pos2=pos1
            ins_seq = str(record.ALT[0])
        else:
            chr2, pos2, ori = get_alt(record)
        qual = record.QUAL
        filt = 'PASS'
        info = modify_info(record, ins_seq)
        form = 'GT'
        genotype = extract_genotype(record)
        
        new_record={'chr1':chr1, 'pos1':pos1, 'ident':ident, 'ref':ref, 'qual':qual, 'filt':filt, 'info':info, 'form':form, 'gt':genotype, 'chr2':chr2, 'pos2':pos2}
        
        
        flag = True #Remove duplicates
        for i, r in new_records_d.items():            
            if (chr1 == r['chr1'] or chr1 == r['chr2']) and (chr2 == r['chr1'] or chr2 == r['chr2']):
                if r['pos1'] == pos2 and r['pos2'] == pos1:
                    flag = False
                
        if flag:
            try:
                if int(chr2) < int(chr1):
                    chr1, chr2 = chr2, chr1
                    pos1, pos2 = pos2,pos1
                elif int(chr2)== int(chr1):
                    if int(pos2) < int(pos1):
                        pos1, pos2 = pos2,pos1
            except ValueError:
                if chr2 == chr1:
                    if int(pos2) < int(pos1):
                        pos1, pos2 = pos2,pos1
            if alt != '<INS>':
                alt = construct_alt(chr2, pos2, ori ,ref)
            new_record_list=[chr1, pos1, ident, ref, alt, qual, filt, info, form, genotype]
            new_record_list = [ str(x) for x in new_record_list]
            new_records_d[counter] = new_record
            counter += 1
            
            print '\t'.join(new_record_list)
        
            
            
            
        


if __name__ == '__main__':
    main()
