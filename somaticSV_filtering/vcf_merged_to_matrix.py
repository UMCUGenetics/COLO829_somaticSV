#! /bin/python
from sys import argv
import vcf

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
    return [pcr, capture,bionano]



with open(argv[1], 'r') as f:
    r = vcf.Reader(f)
    print '\t'.join(["illumina", "pacbio", "tenx", "nanopore","pcr","capture","bionano"])
    for record in r:
        print '\t'.join(get_tech_info(record)+get_validation_info(record))
