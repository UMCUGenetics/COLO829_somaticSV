#! /bin/python
from sys import argv
import vcf

def get_tech_info(record):
    illumina,pacbio,tenx,nanopore = '0'*4
    if 'ILLUMINA' in record.FILTER:
        illumina = '1'
    if 'PACBIO' in record.FILTER:
        pacbio = '1'
    if '10X' in record.FILTER:
        tenx = '1'
    if 'NANOPORE' in record.FILTER:
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
    print '\t'.join(["illumina", "pacbio", "tenx", "nanopore","pcr","capture","bionano","mult_support"])
    for record in r:
	techList = get_tech_info(record)
	valList = get_validation_info(record)
	if sum([int(i) for i in techList]) > 1:
		multSupp = "1"
	else:
		multSupp = "0"
	valList.append(multSupp)
	print '\t'.join(techList + valList)
