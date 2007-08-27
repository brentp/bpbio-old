def split_blast(blastfile,queryfn,subjectfn,qorg,sorg):
    """
    take a single blastfile of m chromosomes vs n chromosomes and split it
    into m * n files based on the queryfn, subjectfn functions which take
    a string accn and return the chromsome: e.g.
    >>> queryfn('At2g26430')
    '2'
    >>> subjectfn('Os01g12345')
    '01'
    """

    outfiles = {}
    for line in open(blastfile):
        line = line.strip().split("\t")
        querychr = queryfn(line[0])
        subjectchr = subjectfn(line[1])
        qchr = str(querychr)
        schr = str(subjectchr)
        if None in (qchr, schr):
            continue
        key = qchr + '_' + schr
        filename = qorg + 'chr' + qchr + '_vs_' \
                 + sorg + 'chr' + schr + '.blast'

        #if not key in outfiles:
        #    outfiles[key] = open(filename,'w')
        #print >>outfiles[key], "\t".join(line)
        for qs in (0, 1):
            if '.' in line[qs]:
                line[qs] = line[qs][:line[qs].find('.')]

        fh = open(filename,'a')
        print >>fh, "\t".join(line)
        fh.close()



if __name__ == "__main__":
    import sys
    import re
    #split_blast(sys.argv[1],lambda x: x[2:3], lambda x: x[1:3], qorg=sys.argv[2], sorg=sys.argv[3])
    def sbchr(accn):
        chr = re.search('SB(\d+)',accn).groups(0)[0]
        if int(chr) > 10:
            return None
        else: 
            return chr
    def oschr(accn):
        return re.search('OS(\d+)',accn).groups(0)[0]

    
    maizefn = lambda x: x.split('|')[1] 
    #split_blast(sys.argv[1], lambda x: x[2:3], lambda x: x[2:3], qorg=sys.argv[2], sorg=sys.argv[3])
    #split_blast(sys.argv[1], qschr, qschr, qorg=sys.argv[2], sorg=sys.argv[3])
    #split_blast(sys.argv[1], lambda x: x.lstrip('>')[2:4], lambda x: x.lstrip('>')[2:4], qorg=sys.argv[2], sorg=sys.argv[3])

    #split_blast(sys.argv[1], maizefn, maizefn, qorg=sys.argv[2], sorg=sys.argv[3])
    split_blast(sys.argv[1], sbchr, sbchr, qorg=sys.argv[2], sorg=sys.argv[3])

 
