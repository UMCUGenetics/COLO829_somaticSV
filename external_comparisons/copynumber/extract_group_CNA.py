#! /bin/python

from sys import argv






def readFile(filename):
    '''Will read the file and return the dictionaries with the CN calls for each group
    Only CN calls with a frequency >= 0.1 are considered'''
    dictA = {}
    dictB = {}
    dictC = {}
    dictD = {}
    with open(filename, 'r') as infile:
        infile.next()
        for line in infile:
            chrom, start, end, cn, event, a, b, c, d = line.rstrip().split('\t')
            if float(a) >= 0.1:
                dictA.setdefault(chrom, [])
                dictA[chrom].append((start, end, cn, a))
            if float(b) >= 0.1:
                dictB.setdefault(chrom, [])
                dictB[chrom].append((start, end, cn, b))
            if float(c) >= 0.1:
                dictC.setdefault(chrom, [])
                dictC[chrom].append((start, end, cn, c))
            if float(d) >= 0.1:
                dictD.setdefault(chrom, [])
                dictD[chrom].append((start, end, cn, d))
            
    return dictA, dictB, dictC, dictD

def remove_overlaps(d):
    '''Removes the CN event overlaps. Takes the one with highest frequency'''
    newD = {}
    for chrom, cnList in d.items():
        print "Chromosome " + chrom
        current_pos = 0
        index_removal=[]
        current_element = cnList[0]
        for index, event in enumerate(cnList):
            start = int(event[0])
            end = int(event[1])
            for event2 in cnList:
                if event == event2:
                    continue
                start2 = int(event2[0])
                end2 = int(event2[1])
                
                if (start > start2 and start < end2) or (start2 > start and start2 < end) or (end > start2 and end < end2) or (end2 > start and end2 < end):
                    if float(event2[3]) > float(event[3]):
                        print "Removing event " + str(event) + "; maintaining event " + str(event2)
                        index_removal.append(index)
                        continue
        for i in sorted(set(index_removal), reverse=True):
            del cnList[i]
        newD[chrom] = cnList
    return newD
                
                
        
def writeFile(filename, d):
    '''Writes filename using d of CN events'''
    with open(filename, 'w') as outfile:
        for chrom, cnList in sorted(d.items()):
            for event in cnList:
                start, end, cn, freq = event
                cn = str(int(cn)-2) #Normalize to 0 for circos plotting
                outfile.write("\t".join(["hs"+chrom, start, end, cn]))
                outfile.write("\n")
            
            
def main():
    all_cna_file = argv[1]
    dictA, dictB, dictC, dictD = readFile(all_cna_file)
    print "Group A"
    dA = remove_overlaps(dictA)
    writeFile("singleCell_groupA_cnv.circos", dA)
    print "Group B"
    dB = remove_overlaps(dictB)
    writeFile("singleCell_groupB_cnv.circos", dB)

    print "Group C"
    dC = remove_overlaps(dictC)
    writeFile("singleCell_groupC_cnv.circos", dC)

    print "Group D"
    dD = remove_overlaps(dictD)
    writeFile("singleCell_groupD_cnv.circos", dD)

    
    
    
        




if __name__ == "__main__":
    main()