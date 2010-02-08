import sys
import math
from pyfasta import Fasta
import numpy as np

hists = (
    "GSM310840_Cy5_Cy3.txt",
    "GSM310841_Cy5_Cy3.txt",
    "GSM310842_Cy5_Cy3.txt",
    "GSM310843_Cy5_Cy3.txt",
)

posns = (
    "GPL7143_2007-06-19_ATH1_chr_all_meth01.pos",
    "GPL7143_2007-06-19_ATH1_chr_all_meth02.pos",
    "GPL7143_2007-06-19_ATH1_chr_all_meth03.pos",
)


"""
SEQ_ID  CHROMOSOME  PROBE_ID    POSITION    LENGTH  COUNT
CHR1v01212004   chr1v01212004   CHR1V01212004FS008761566    8761566 46  1
CHR1v01212004   chr1v01212004   CHR1V01212004FS013192017    13192017    52  1
"""
def parse_posns(pfiles):

    lookups = {}
    for f in pfiles:
        header = None
        for line in open(f):
            line = line.rstrip().split("\t")
            if header is None:
                header = line 
                continue
            lookups[line[2]] = int(line[3]), int(line[4])
    return lookups


"""
ID_REF  Cy5 Cy3
CHR1V01212004FS000000001    39524.11    31755.55
CHR1V01212004FS000000061    3463.56 3682.89
CHR1V01212004FS000000121    1434.22 1427.44
"""
import collections
def parse_hists(fhists):

    vals_by_id = collections.defaultdict(list)
    for f in fhists:
        header = None
        for line in open(f):
            line = line.rstrip().split("\t")
            if header is None:
                header = line 
                continue
            id, cy5, cy3 = line[0], float(line[1]), float(line[2])
            vals_by_id[id].append((cy5, cy3))

    # TODO: see zscore here: for confidence in measurement...
    # http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM310842
    def avg(li): return sum(li)/len(li)
    for id, c53li in vals_by_id.iteritems():
        cy5 = avg([x[0] for x in c53li])
        cy3 = avg([x[1] for x in c53li])
        yield id, math.log(1 + (cy5 / cy3))


if __name__ == "__main__":
    plookups = parse_posns(posns)
    fasta = Fasta("/opt/src/flatfeature/data/thaliana_v9.fasta")

    arrs = {}
    counts = {}
    for seqid, seq in fasta.iteritems():
        arrs[seqid] = np.zeros((len(seq),), dtype=np.float32)
        counts[seqid] = np.zeros((len(seq),), dtype=np.float32)
    print arrs.keys()

    success = 0
    for id, val in parse_hists(hists):
        start, length = plookups[id]
        #print id, val, start, length
        seqid = id[3]
        if not seqid in "12345": continue
        a = arrs[seqid][start - 1: start + length]
        counts[seqid][start - 1: start + length] += 1
        c = counts[seqid][start -1: start + length]
    

        # keep the average. since there are overlaps, this 
        # weights by the number of existing measurements taht have
        # already contributed to the value.
        arrs[seqid][start - 1: start + length] = ((c - 1) * a + val) / c
        success += 1
    print "ADDED:", success
    for seqid, arr in arrs.iteritems():
        arr.tofile("histone.%s.bin" % seqid)
