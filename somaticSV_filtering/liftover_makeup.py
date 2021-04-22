from sys import argv
import re

#Truthset was liftover from GRCh37 to GRCh38 using ENSEMBL Assembly Converter
#https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter
#Problem is that liftover does not change the ALT field, so that's what we do here. 
#I have not checked orientation though


#store the breakpoints in pairs
d = {}
with open(argv[1], 'r') as invcf:
	for line in invcf:
		line = line.rstrip()
		if line.startswith("#"):
			print line
			continue
		chrm, pos1, fullid, ref, alt_orig, qual, filt, info, gt1, gt2 = line.split("\t")
		k = "_".join(fullid.split("_")[0:2])
		if k in d:
			d[k][2] = {"chrm":chrm, "pos1":pos1,"fullid":fullid, "ref":ref, "alt_orig":alt_orig,
					"qual":qual, "filt":filt, "info":info, "gt1":gt1, "gt2":gt2}
		else:
			d[k] = {1:{"chrm":chrm, "pos1":pos1,"fullid":fullid, "ref":ref, "alt_orig":alt_orig,
					"qual":qual, "filt":filt, "info":info, "gt1":gt1, "gt2":gt2}}


#Change the second coordinate and the SVLEN info field
for k, bp_d in d.items():
	for i in [1,2]:
		if i == 1:
			i2 = 2
		else:
			i2 = 1

		alt1 = bp_d[i]["alt_orig"]
		info1 = bp_d[i]["info"]

		if  alt1 == "<INS>":
			bp_d[i]["alt_new"] = bp_d[i]["alt_orig"]
			bp_d[i]["info_new"] = bp_d[i]["info"]
			break

		alt_match = re.search("^(\w*\]|\w*\[)(\w+):(\d+)(\]\w*|\[\w*)$", alt1)
		ori1 = alt_match.group(1)
		chrm_orig = alt_match.group(2)
		pos_orig = alt_match.group(3)
		ori2 = alt_match.group(4)

		chrm_new = bp_d[i2]["chrm"]
		pos_new = bp_d[i2]["pos1"]
		alt_new = ori1+chrm_new+":"+pos_new+ori2

		#Change SVLEN since it changes for a couple of SVs
		if chrm_new == bp_d[i]["chrm"]:
			new_svlen = str(abs(int(pos_new) - int(bp_d[i]["pos1"])))
		else: 
			new_svlen = "0"
		info_new = re.sub('SVLEN=[0-9]+;', bp_d[i]["info"], "SVLEN="+new_svlen+";")

		bp_d[i]["alt_new"] = alt_new
		bp_d[i]["info_new"] = info_new

for bp, bp_d in d.items():
	for i in [1,2]:
		try:
			print "\t".join([bp_d[i]["chrm"],bp_d[i]["pos1"],bp_d[i]["fullid"],bp_d[i]["ref"],
			bp_d[i]["alt_new"],bp_d[i]["qual"],bp_d[i]["filt"],bp_d[i]["info_new"],bp_d[i]["gt1"],bp_d[i]["gt2"]])
		except KeyError:	#For the INS breakpoints there is no second one
			continue








