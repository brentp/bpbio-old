

import sys


def add_locs_to_blast():
    blastfile = sys.argv[1]
    sep = len(sys.argv) > 2 and sys.argv[2] or '|'

    for line in open(blastfile):
        line = line.strip().split('\t')
        q = line[0].split(sep) # name|start|stop
        s = line[1].split(sep) # name|start|stop
        locs = map(int,line[6:10])
        locs[0]+=int(q[1])
        locs[1]+=int(q[1])
        locs[2]+=int(s[1])
        locs[3]+=int(s[1])
        line[6:10] = map(str,locs)
        print "\t".join(line)
