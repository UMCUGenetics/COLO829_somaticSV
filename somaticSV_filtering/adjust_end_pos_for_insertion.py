#!/usr/bin/python

import sys
import re

with(open(sys.argv[1],'r')) as vcf:
    for line in vcf:
	line = line.rstrip()
	if line.startswith("#"):
	    print(line)
	else:
	    columns = line.split("\t")
	    end_match = re.search(";END=(\d+);",columns[7])
	    end = None
	    if end_match:
	        end = end_match.group(1)
	    if columns[1] == end:
	        end_adj = int(end)+1
	        columns[7] = columns[7].replace(";END="+end,";END="+str(end_adj))
	        print("\t".join(columns))
	    else:
	        print(line)
