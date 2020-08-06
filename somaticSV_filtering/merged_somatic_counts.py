#! /bin/python
import vcf as pyvcf

def get_sv_info( record ):
    chr1 = str(record.CHROM)
    pos1 = str(record.POS)
    svtype = str(record.INFO['SVTYPE_INF'])
    svlen = str(record.INFO['SVLEN'])
    ori = str(record.INFO['STRANDS'])
    if svtype == "INS":
        chr2 = chr1
        pos2 = str(int(pos1)+int(record.INFO['SVLEN']))
    else:
        chr2 = str(record.ALT[0].chr)
    	pos2 = str(record.ALT[0].pos)
    sv_info = [str(record.ID),chr1,pos1,chr2,pos2,svtype,ori,svlen]
    return(sv_info)

def get_tech_info(record):
    illumina,pacbio,tenx,nanopore = '0'*4
    if 'illumina' in record.INFO['SUPP_SEQ']:
        illumina = '1'
    if 'pacbio' in record.INFO['SUPP_SEQ']:
        pacbio = '1'
    if '10X' in record.INFO['SUPP_SEQ']:
        tenx = '1'
    if 'nanopore' in record.INFO['SUPP_SEQ']:
        nanopore = '1'
    return [illumina, pacbio, tenx, nanopore]
    
def get_validation_info(record):
    pcr,capture,bionano = '0'*3
    if 'PCR' in record.FILTER:
        pcr = "1"
    if 'CAPTURE' in record.FILTER:
        capture = "1"
    if 'BIONANO' in record.FILTER:
        bionano = "1"
    return [pcr, capture, bionano]


def main():
    mult = "somatic.multiplesupport.validations.vcf"
    one = "somatic.onesupport.validations.vcf"
    multout="somatic.multiplesupport.validations.txt"
    oneout = "somatic.onesupport.validations.txt"

    mult_reader = pyvcf.Reader(open(mult, 'r'))    
    with open(multout, 'w') as outf:
        outf.write("\t".join([
    "svid","chr1","pos1","chr2","pos2","svtype","ori","svlen",
    "illumina","pacbio","tenx","nanopore","pcr", "capture", "bionano"])+'\n')
        for record in mult_reader:
            sv_info = get_sv_info(record)
            tech_info = get_tech_info(record)
            val_info = get_validation_info(record)
            outf.write("\t".join([
                "\t".join(sv_info),
                "\t".join(tech_info),
                "\t".join(val_info)]) + "\n")
        
        
    one_reader = pyvcf.Reader(open(one, 'r'))    
    with open(oneout, 'w') as outf:
        outf.write("\t".join([
    "svid","chr1","pos1","chr2","pos2","svtype","ori","svlen",
    "illumina","pacbio","tenx","nanopore","pcr", "capture", "bionano"])+'\n')
        for record in one_reader:
            sv_info = get_sv_info(record)
            tech_info = get_tech_info(record)
            val_info = get_validation_info(record)
            outf.write("\t".join([
                "\t".join(sv_info),
                "\t".join(tech_info),
                "\t".join(val_info)]) + "\n")


if __name__ == "__main__":
    main()



