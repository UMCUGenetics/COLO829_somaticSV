#! /usr/bin/python

import argparse
import re

file_types = dict()
file_types["vcf"] = [ "delly", "manta", "pindel", "sniffles", "lumpy", "nanosv", "haplotypecaller", "1kgsv", "1kgindel" ]
file_types["bed"] = [ "all types of bed files"]
file_types["other"] = [ "freec", "mobster" ]

structural_variants = dict()
bed = dict()
size_cutoff = 1000
reciprocal_overlap = 0.7

parser = argparse.ArgumentParser(description='Annotate a sv vcf file with another file.')
parser.add_argument('--input', nargs='?', help='vcf file to annotate with <file2>', required=True)
parser.add_argument('--file2', nargs='?', help='use this file for annotation', required=True)
parser.add_argument('--distance', nargs='?', help='set distance to extend the confidence interval.')
parser.add_argument('--annotation', nargs='?', help='use this name for annotation in stead of "intersect".', default="intersect")
parser.add_argument('--both_sides', action='store_true', default=False, help='report intersect if both sides of the sv are overlapping the feature in bed file.')
parser.add_argument('--one_feature', action='store_true', default=False, help='report intersect if start and end of the feature in bed file overlaps with confidence interval of one sv.')
args = parser.parse_args()

def parse_bed_file( bed_file ):
        i = 1
        with open(bed_file) as b:
                for line in b:
                        line = line.rstrip()
                        if not line.startswith("#"):
                                columns = line.split("\t")
                                chrom, start, end = columns[:3]
                                if len(columns) >= 4:
                                        i = columns[3]
                                chrom = chrom.replace("chr","")
                                start = int(start)
                                end = int(end)
                                if not chrom in bed:
                                        bed[chrom] = dict()
                                if not start in bed[chrom]:
                                        bed[chrom][start] = dict()
                                        
                                if not ( chrom or start or end ):
                                        sys.exit("Unknown bed file format")
                                bed[chrom][start][end] = str(i)
                                i = int(i)+1

def parse_other_file ( file ):
        with open(file) as f:
                for line in f:
                        line = line.rstrip()
                        if not line.startswith("#") and not line.startswith("Chr\tMobile Element"):
                                columns = line.split("\t")
                                chrom, pos1, pos2, svtype, id = 0, 0, 0, "", ""
                                pos1a, pos1b = pos1, pos1
                                pos2a, pos2b = pos2, pos2
                                
                                # parse input
                                ## freec
                                if len(columns) == 5 and ( columns[4] == "gain" or columns[4] == "loss" ):
                                        chrom, pos1, pos2 = columns[:3]
                                        pos1 = int(pos1)
                                        pos2 = int(pos2)
                                        chrom2 = chrom
                                        svtype = columns[4]
                                        pos1a, pos1b = pos1, pos1
                                        pos2a, pos2b = pos2, pos2
                                        pos1a -= 1000
                                        pos2b += 1000
                                        id = "_".join([columns[4],chrom,str(pos1),str(pos2)])
                                ## mobster
                                if len(columns) == 28 and ( columns[1] == "ALU" or columns[1] == "L1" or columns[1] == "SVA" ):
                                        chrom = columns[0]
                                        pos1 = int(columns[2])
                                        pos2 = int(columns[2])
                                        svtype = "INS"
                                        chrom2 = chrom
                                        pos1a = int(columns[3])
                                        pos1b = pos1
                                        pos2a = pos2
                                        pos2b = int(columns[4])
                                        id = "_".join([columns[1],chrom,str(pos1),str(pos2)])
                                
                                chrom = chrom.replace("chr","")
                                chrom2 = chrom2.replace("chr","")
                                
                                # set svtype
                                ## freec
                                if svtype == "loss":
                                        svtype = "DEL"
                                elif svtype == "gain":
                                        svtype = "DUP"
                                
                                # set orientation
                                orientation = ""
                                if svtype == "DEL":
                                        orientation = "TH"
                                if svtype == "INS":
                                        orientation = "TH"
                                elif svtype == "DUP":
                                        orientation = "HT"
                                
                                # add second orientation to structural variant
                                orientations = [ orientation ]
                                if svtype == "INS":
                                        orientation2 = "HT"
                                        orientations.append( orientation2 )
                                
                                if not ( chrom or chrom2 or pos1 or pos2 or orientations ):
                                        sys.exit("Unknown other type of file format")
                                
                                for ori in orientations:
                                        # switch coordinates from lowest to highest
                                        if pos1 > pos2 and chrom == chrom2:
                                                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
                                                ori = ori[::-1]
                                        elif chrom2 < chrom:
                                                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
                                                ori = ori[::-1]
                                        
                                        sv_key = "\t".join([chrom, chrom2])
                                        info_value = "\t".join([str(pos1), str(pos2), svtype, id])
                                        
                                        if not sv_key in structural_variants:
                                                structural_variants[sv_key] = dict()
                                        if not pos1a in structural_variants[sv_key]:
                                                structural_variants[sv_key][pos1a] = dict()
                                        if not pos1b in structural_variants[sv_key][pos1a]:
                                                structural_variants[sv_key][pos1a][pos1b] = dict()
                                        if not pos2a in structural_variants[sv_key][pos1a][pos1b]:
                                                structural_variants[sv_key][pos1a][pos1b][pos2a] = dict()
                                        structural_variants[sv_key][pos1a][pos1b][pos2a][pos2b] = info_value

def parse_bpi_tsv (file):
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("seqnames1"): 
                continue
            columns = line.rstrip().split('\t')            
            chrom, pos1, end1, strand1, chrom2, pos2, end2, strand2, score, svtype, bafup, bafdown, depth, sample = columns
            chrom = chrom.replace("chr","")
            chrom2 = chrom2.replace("chr","")
            pos1 = int(pos1)
            pos2 = int(pos2)
            orientation = ""
            if svtype == "DEL":
                orientation = "TH"
            elif svtype == "INS":
                orientation = "TH"
            elif svtype == "DUP":
                orientation = "HT"
            elif svtype == "BND":
                if strand1 == "-":
                    orientation = "H"
                elif strand1 == "+":
                    orientation = "T"
                if strand2 == "-":
                    orientation += "H"
                elif strand2 == "+":
                    orientation += "T"
            elif svtype.startswith("INV"):
                orientation = svtype.replace("(", "").replace(")", "").replace("2", "").upper()
                svtype = "INV"
            
            # set confidence interval around start
            pos1a, pos1b = pos1, pos1
            pos2a, pos2b = pos2, pos2
            if pos1 > pos2 and chrom == chrom2:
                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
            elif chrom2 < chrom:
                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
                
                
            sv_key = "\t".join([chrom, chrom2])
            ident = svtype+score
            info_value = "\t".join([str(pos1), str(pos2), svtype, ident])  
            if not sv_key in structural_variants:
                structural_variants[sv_key] = dict()
            if not pos1a in structural_variants[sv_key]:
                structural_variants[sv_key][pos1a] = dict()
            if not pos1b in structural_variants[sv_key][pos1a]:
                structural_variants[sv_key][pos1a][pos1b] = dict()
            if not pos2a in structural_variants[sv_key][pos1a][pos1b]:
                structural_variants[sv_key][pos1a][pos1b][pos2a] = dict()
            structural_variants[sv_key][pos1a][pos1b][pos2a][pos2b] = info_value
    #print structural_variants
    #p1, p2, t, i = structural_variants[sv_key][p1a][p1b][p2a][p2b].split("\t")



def parse_vcf_file ( vcf_file ):
        with open(vcf_file) as vcf:
                for line in vcf:
                        line = line.rstrip()
                        if line.startswith("##"):
                                if vcf_file == args.input:
                                        print( line )
                        elif line.startswith("#"):
                                if vcf_file == args.input:
                                        print( "##FILTER=<ID="+args.annotation+",Description=\"Overlap with "+args.file2+" within distance of "+args.distance+"\">" )
                                        print( line )
                        else:
                                intersect = False
                                intersect_ids = list()
                                
                                # parse vcf fields
                                columns = line.split("\t")
                                chrom, pos1, id, ref, alt, qual, filter, info = columns[:8]
                                chrom = chrom.replace("chr","")
                                pos1 = int(pos1)
                                
                                # create regex match varianbles
                                alt_match = re.search("^(\w*\]|\w*\[)(\w+):(\d+)(\]\w*|\[\w*)$", alt)
                                svtype_match = re.search("SVTYPE=(\w+)(;|$)",info)
                                chrom2_match = re.search("CHR2=(.+?)(;|$)", info)
                                end_match = re.search("END=(\d+)(;|$)", info)
                                ct_match = re.search("CT=(\d+)to(\d+)(;|$)",info) # DELLY
                                strands_match = re.search("STRANDS=([\+-])([\+-])(;|$)", info) # SNIFFLES AND LUMPY
                                
                                # set svtype
                                svtype = ""
                                if svtype_match:
                                        svtype = svtype_match.group(1)
                                ## 1kgsv specific notation
                                if svtype == "HERV" or svtype == "SVA" or svtype[0:4] == "LINE" or svtype[0:3] == "ALU":
                                        svtype = "INS"
                                
                                ## 1kgindel and haplotypecaller specific notation
                                elif "INDEL" in info:
                                        if (len(alt)-len(ref))+1 > 0:
                                                svtype = "INS"
                                        elif (len(alt)-len(ref))+1 < 0:
                                                svtype = "DEL"
                                #1kgsv specific notation
                                elif svtype == "CNV":
                                        if "DUP" in info:
                                                svtype = "DUP"
                                        elif "DEL" in info:
                                                svtype = "DEL"
                                ##haplotypecaller
                                #else:
                                        #if (len(alt)-len(ref))+1 > 0:
                                                #svtype = "INS"
                                        #elif (len(alt)-len(ref))+1 < 0:
                                                #svtype = "DEL"
                        
                                # set second chromosome
                                chrom2 = chrom
                                if chrom2_match:
                                        chrom2 = chrom2_match.group(1)
                                elif alt_match:
                                        chrom2 = alt_match.group(2)
                                chrom2 = chrom2.replace("chr","")
                                
                                # set end position
                                pos2 = 0
                                if end_match:
                                        pos2 = end_match.group(1)
                                elif alt_match:
                                        pos2 = alt_match.group(3)
                                ## 1kgindel and haplotypecaller specific notation
                                elif "INDEL" in info:
                                        if ((len(alt)-len(ref))+1) < 0:
                                                pos2 = pos1-(len(alt)-len(ref))+1
                                        else:
                                                pos2 = pos1
                                elif svtype == "INS":
                                        pos2 = pos1+1
                                pos2 = int(pos2)
                                
                                # set orientation
                                orientation = ""
                                ## set orientation for interchromosomal based on ALT field
                                if alt_match:
                                        if alt_match.group(1) == "]":
                                                orientation = "HT"
                                        elif alt_match.group(1) == "[":
                                                orientation = "HH"
                                        elif alt_match.group(4) == "]":
                                                orientation = "TT"
                                        elif alt_match.group(4) == "[":
                                                orientation = "TH"

                                ## set orientation for intrachromosomal sv's
                                elif svtype == "DEL":
                                        orientation = "TH"
                                elif svtype == "INS":
                                        orientation = "TH"
                                elif svtype == "DUP":
                                        orientation = "HT"
                                ## delly specific notation
                                elif ct_match: 
                                        for g in ct_match.groups()[:2]:
                                                if int(g) == 5:
                                                        orientation += "H"
                                                elif int(g) == 3:
                                                        orientation += "T"
                                ## manta specific notation
                                elif "INV5" in info:
                                        orientation = "HH"
                                elif "INV3" in info:
                                        orientation = "TT"
                                # pindel specific notation
                                elif svtype == "RPL":
                                        orientation = "TH"
                                # sniffles and lumpy specific notation
                                elif strands_match:
                                        for g in strands_match.groups()[:2]:
                                                if g == "-":
                                                        orientation += "H"
                                                elif g == "+":
                                                        orientation += "T"
                                
                                # set confidence interval around start
                                pos1a, pos1b = pos1, pos1
                                
                                if "CIPOS" in info:
                                        cipos = map(int, re.search("CIPOS=-*(\d+,\d+)", info).group(1).split(","))
                                        pos1a -= cipos[0]
                                        pos1b += cipos[1]
                                
                                # set confidence interval around end
                                pos2a, pos2b = pos2, pos2
                                if "CIEND" in info:
                                        ciend = map(int, re.search("CIEND=-*(\d+,\d+)", info).group(1).split(","))
                                        pos2a -= ciend[0]
                                        pos2b += ciend[1]
                                
                                # add second orientation to structural variant
                                orientations = [ orientation ]
                                if svtype == "INS":
                                        orientation2 = "HT"
                                        orientations.append( orientation2 )
                                ## pindel, lumpy and 1kgsv specific INV
                                elif svtype == "INV" and not orientation:
                                        orientation = "TT"
                                        orientations.append( orientation )
                                        orientation2 = "HH"
                                        orientations.append( orientation2 )
                                ## sniffles specific DUP/INS
                                elif svtype == "DUP/INS":
                                        orientation = "HT"
                                        orientations.append( orientation )
                                        orientation2 = "TH"
                                        orientations.append( orientation2 )
                                ## sniffles specific DEL/INV
                                elif svtype == "DEL/INV":
                                        orientation = "TH"
                                        orientations.append( orientation )
                                        orientation2 = "TT"
                                        orientations.append( orientation2 )
                                        orientation3 = "HH"
                                        orientations.append( orientation3 )
                                        
                                # rename INV and TRA to BND
                                if svtype == "INV" or svtype == "TRA":
                                        svtype = "BND"
                                
                                # add extra distance to the confidence interval
                                if args.distance and vcf_file == args.input:
                                        pos1a -= int(args.distance)
                                        pos1b += int(args.distance)
                                        pos2a -= int(args.distance)
                                        pos2b += int(args.distance)					
                                
                                if not ( chrom or chrom2 or pos1 or pos2 or orientations ):
                                        sys.exit("Unknown vcf file format")
                                
                                for ori in orientations:
                                        # switch coordinates from lowest to highest
                                        if pos1 > pos2 and chrom == chrom2:
                                                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
                                                ori = ori[::-1]
                                        elif chrom2 < chrom:
                                                chrom, pos1, pos1a, pos1b, chrom2, pos2, pos2a, pos2b = chrom2, pos2, pos2a, pos2b, chrom, pos1, pos1a, pos1b
                                                ori = ori[::-1]
                                        
                                        sv_key = "\t".join([chrom, chrom2])
                                        info_value = "\t".join([str(pos1), str(pos2), svtype, id])
                                        
                                        # if file is the annotation file
                                        if vcf_file == args.file2:
                                                if not sv_key in structural_variants:
                                                        structural_variants[sv_key] = dict()
                                                if not pos1a in structural_variants[sv_key]:
                                                        structural_variants[sv_key][pos1a] = dict()
                                                if not pos1b in structural_variants[sv_key][pos1a]:
                                                        structural_variants[sv_key][pos1a][pos1b] = dict()
                                                if not pos2a in structural_variants[sv_key][pos1a][pos1b]:
                                                        structural_variants[sv_key][pos1a][pos1b][pos2a] = dict()
                                                if not pos2b in structural_variants[sv_key][pos1a][pos1b][pos2a]:
                                                    structural_variants[sv_key][pos1a][pos1b][pos2a][pos2b] = list()
                                                structural_variants[sv_key][pos1a][pos1b][pos2a][pos2b].append(info_value)
                                        
                                        # if file is the input file
                                        elif vcf_file == args.input:
                                                if structural_variants and sv_key in structural_variants:                                                        
                                                        for p1a in structural_variants[sv_key]:
                                                                for p1b in structural_variants[sv_key][p1a]:
                                                                        #print "CI 1 CHECK"
                                                                        #print "\t".join([str(p1a), str(p1b), str(pos1), str(pos1a), str(pos1b)])
                                                                        # next if confidence intervals on the 'right' side do not overlap
                                                                        if not (p1a <= pos1b and p1b >= pos1a):
                                                                                continue
                                                                        #print "PASS CHECK 1"
                                                                        for p2a in structural_variants[sv_key][p1a][p1b]:
                                                                                for p2b in structural_variants[sv_key][p1a][p1b][p2a]:
                                                                                    for l in structural_variants[sv_key][p1a][p1b][p2a][p2b]:
                                                                                                # next if confidence intervals on the 'left' side do not overlap
                                                                                                if not (p2a <= pos2b and p2b >= pos2a):
                                                                                                        continue
                                                                                                #print "PASS CHECK 2"
#												p1, p2, t, i = structural_variants[sv_key][p1a][p1b][p2a][p2b].split("\t")
                                                                                                p1, p2, t, i = l.split("\t")
                                                                                                #if (i == id):
                                                                                                    #continue
                                                                                                #print "PASS CHECK 3"
                                                                                                p1 = int(p1)
                                                                                                p2 = int(p2)
                                                                                                # Sv's must overlap if it is < <size_cutoff> and one of the sv's is not an insertion
                                                                                                if ( chrom == chrom2 and (pos2-pos1)<size_cutoff and ( svtype != "INS" and t != "INS" ) ):
                                                                                                        if ( p1 <= pos2 and p2 >= pos1 ):
                                                                                                                
                                                                                                                # insertion and duplication only have to overlap
                                                                                                                if ( ( t == "INS" and svtype == "DUP" ) or ( t == "DUP" and svtype == "INS" ) ):
                                                                                                                        intersect = True
                                                                                                                        intersect_ids.append(i)
                                                                                                                # other sv type needs a minimal reciprocal overlap
                                                                                                                else:
                                                                                                                        len1 = pos2-pos1
                                                                                                                        len2 = p2-p1
                                                                                                                        if len1 == 0:
                                                                                                                                len1 = 1
                                                                                                                        if len2 == 0:
                                                                                                                                len2 = 1
                                                                                                                        sorted_coords = [pos1, pos2, p1, p2]
                                                                                                                        sorted_coords.sort()
                                                                                                                        olen = sorted_coords[2]-sorted_coords[1]
                                                                                                                        if ( olen/float(len1) > reciprocal_overlap and olen/float(len2) > reciprocal_overlap ):
                                                                                                                                intersect = True
                                                                                                                                intersect_ids.append(i)
                                                                                                else:
                                                                                                        intersect = True
                                                                                                        intersect_ids.append(i)
                                                elif bed:
                                                        if args.one_feature:
                                                                if chrom in bed and chrom == chrom2:
                                                                        for p1 in bed[chrom]:
                                                                                if not ( p1 <= pos1b and p1 >= pos1a ):
                                                                                        continue
                                                                                for p2 in bed[chrom][p1]:
                                                                                        if not ( p2 <= pos2b and p2 >= pos2a ):
                                                                                                continue
                                                                                        intersect = True
                                                                                        intersect_ids.append(bed[chrom][p1][p2])
                                                                                        
                                                        elif args.both_sides:
                                                                intersect1 = False
                                                                intersect2 = False
                                                                intersect_ids1 = list()
                                                                intersect_ids2 = list()
                                                                if chrom in bed:
                                                                        for p1a in bed[chrom]:
                                                                                for p1b in bed[chrom][p1a]:
                                                                                        if not (p1a <= pos1b and p1b >= pos1a):
                                                                                                continue
                                                                                        intersect1 = True
                                                                                        intersect_ids1.append(bed[chrom][p1a][p1b])
                                                                if chrom2 in bed:
                                                                        for p2a in bed[chrom2]:
                                                                                for p2b in bed[chrom2][p2a]:
                                                                                        if not (p2a <= pos2b and p2b >= pos2a):
                                                                                                continue
                                                                                        intersect2 = True
                                                                                        intersect_ids2.append(bed[chrom2][p2a][p2b])
                                                                if intersect1 and intersect2:
                                                                        intersect = True
                                                                        intersect_ids = list(set().union(intersect_ids1,intersect_ids2))
                                                        else:
                                                                if chrom in bed:
                                                                        for p1a in bed[chrom]:
                                                                                for p1b in bed[chrom][p1a]:
                                                                                        if not (p1a <= pos1b and p1b >= pos1a):
                                                                                                continue
                                                                                        intersect = True
                                                                                        intersect_ids.append(bed[chrom][p1a][p1b])
                                                                if chrom2 in bed:
                                                                        for p2a in bed[chrom2]:
                                                                                for p2b in bed[chrom2][p2a]:
                                                                                        if not (p2a <= pos2b and p2b >= pos2a):
                                                                                                continue
                                                                                        intersect = True
                                                                                        intersect_ids.append(bed[chrom2][p2a][p2b])																
                                                                                        
                                
                                if vcf_file == args.input:
                                        if intersect:
                                                if columns[6] == "PASS":
                                                        columns[6] = args.annotation
                                                else:
                                                        columns[6] += ";"+args.annotation
                                                
                                                columns[7] += ";"+args.annotation+"_IDS="+",".join(intersect_ids)
                                        print( "\t".join(columns))
                                        
if ".vcf" in args.file2:
        parse_vcf_file( args.file2 )
elif ".bed" in args.file2:
        parse_bed_file( args.file2 )
elif "bpi.tsv" in args.file2:
        parse_bpi_tsv( args.file2 )
else:
        parse_other_file( args.file2)
parse_vcf_file( args.input )

