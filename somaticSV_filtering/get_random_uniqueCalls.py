#! /bin/python
import random

input="somatic.onesupport.validations.vcf"

d={"nanopore":[],"10X":[],"illumina":[],"pacbio":[]}
illumina_ids = []
pacbio_ids = []
with open(input,'r') as invcf:
    for l in invcf:
        l=l.rstrip()
        if l.startswith("#"):
            print(l)
            continue
        lsplit=l.split('\t')
	if lsplit[4] == "<INS>":
            continue
        if lsplit[6] != "PASS": #Validated
            continue
        tech = lsplit[7].split(";")[1].split("=")[1]
        if tech == "illumina":
            ident = lsplit[2][:-1] #Remove last letter
            if ident in illumina_ids:
                continue
            else:
                illumina_ids.append(ident)
        elif tech == "pacbio":
            ident = lsplit[2]
            if ident in pacbio_ids:
                continue
            else:
                pacbio_ids.append(ident)            
        d[tech].append(l)
        
selected = random.sample(d["nanopore"],80)
for t in ["pacbio", "10X", "illumina"]:
    selected = selected + random.sample(d[t], 40)
for l in selected:
    print l
 
        
