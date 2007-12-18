"""
05||003790001||003800000||10KMER    02||031690001||031700000||10KMER 84.62   78  12  0   8631    8708    7448    7371    3e-07   60.0
05||003790001||003800000||10KMER    02||031690001||031700000||10KMER 84.42   77  12  0   8969    9045    7447    7371    1e-06   58.0
05||003790001||003800000||10KMER    02||026520001||026530000||10KMER 94.81   77  4   0   8993    9069    4091    4015    1e-25    121
05||003790001||003800000||10KMER    02||026520001||026530000||10KMER 88.50   113 11  2   8657    8768    4090    3979    6e-21    105
"""

def updatefile(blastfile, outfile=None):
    if out is None: outfile = blastfile.replace(".blast",".locs.blast")
    fh = open(blastfile)
    out = open(outfile,'w')
    for line in fh:
        line = line.split("\t")
        qchrstart = int(line[0].split("||")[1])
        schrstart = int(line[1].split("||")[1])
        locs = map(int, line[6:10])

        locs[0] += qchrstart # qstart
        locs[1] += qchrstart # qstop
        locs[2] += schrstart # sstart
        locs[3] += schrstart # sstop

        line[6:10] = map(str, locs)
        print >>out, "\t".join(line),
    out.close()
    fh.close()
    return outfile


if __name__ == "__main__":
    import sys
    print updatefile(sys.argv[1])
