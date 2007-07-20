"""
take a single blastfile of m chromosomes vs n chromosomes and split it
into m * n files based on the queryfn, subjectfn functions which take
a string accn and return the chromsome: e.g.
"""
def queryfn(accn):
    """>>> queryfn('At2g26430')
    '2'
    """
    return accn[2]
        

def split_blast(blastfile,queryfn,subjectfn):
    outfiles = {}
    for line in open(blastfile):
        line = line.strip().split("\t")
        querychr = queryfn(line[0])
        subjectchr = subjectfn(line[1])
        key = str(querychr) + '_' + str(subjectchr)
        filename = blastfile + '_' + key
        if not key in outfiles:
            outfiles[key] = open(filename,'w')
        print >>outfiles[key], "\t".join(line)



if __name__ == "__main__":
    import sys
    split_blast(sys.argv[1],lambda x: x[2:4], lambda x: 11)
