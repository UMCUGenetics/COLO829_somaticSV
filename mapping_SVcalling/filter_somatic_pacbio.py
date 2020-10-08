from sys import argv

with open(argv[1], 'r') as infile:
	for line in infile:
		line = line.rstrip()
		if line.startswith("#"):
			print line
		else:
			ref = line.split('\t')[9]
			gt = ref.split(':')[0]
			som_flag = False
			if gt == "./." or gt == "0/0":
				som_flag = True
			info = line.split('\t')[7]
			l_flag = False
			if "BND" in info: 
				l_flag=True
			else:
				svlen=0
				for field in info.split(';'):
					if "SVLEN" in field:
						svlen=abs(int(field.split('=')[-1]))
						break
				if svlen >= 30:
					l_flag = True
			filter = line.split('\t')[6]
			f_flag = False
			if filter == "PASS": #Remove Decoy breakpoints
				f_flag=True
			if som_flag and l_flag and f_flag:
				print line
						
			

