from sys import argv

aln = argv[1]
bin_length = 1000
tech, sample=aln.split('.')[0:2]


with open(aln, 'r') as infile:
    for line in infile:
        if not line.startswith('0'):
            continue
        line = line.rstrip().split('\t')
        ident = line[1]
        l = int(line[2].split(".")[0])
        ident = line[1]
        print ','.join([tech, sample, ident, str(l)])
                    
