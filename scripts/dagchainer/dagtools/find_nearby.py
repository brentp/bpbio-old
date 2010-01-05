"""
the diags_file contains the hits that dagchainer found, the all_file is the full list of blast hits.
this script goes through all of the hits in the dag file and adds any hit from the all_file that is within 'dist' of 
any hit in the diag.
"""
from scipy.spatial import cKDTree
import numpy as np
import sys
from cdagline import DagLine
import collections


def read_dag_to_tree(all_hits):
    """create an rtree, using query as x, subject as y
    do this for all htis, then for each diag, do an intersection
    (+bbuffer) to find nearby 
    """
    trees = collections.defaultdict(list)
    lines = collections.defaultdict(list)
    for sline in open(all_hits):
        if sline[0] == '#': continue
        line = sline[:-1].split("\t")
        chrs = (line[0], line[4])
        # so save the index, which will return i when queried
        # and associate i with the text line vie the dict.
        trees[chrs].append((int(line[2]),  int(line[6])))
        lines[chrs].append(line)

    for chrs in trees:
        trees[chrs] = cKDTree(np.array(trees[chrs]))
    return trees, lines
        
def main(dist, diags, all_hits):
    """empty docstring"""
    trees, lines = read_dag_to_tree(all_hits)
    seen = {}
    for sline in open(diags):
        if sline[0] == '#': 
            # reset seen for each new diagonal...
            # seen = {}
            print sline.rstrip()
            continue
        d = DagLine(sline)
        chrs = (d.a_seqid, d.b_seqid)
        seen[(d.a_accn, d.b_accn)] = True

        tree = trees[chrs]
        info_lines = lines[chrs]
        
        q, s = d.a_start, d.b_start

        dists, idxs = tree.query([q, s], p=2, distance_upper_bound=dist, k=64)
        # values set to tree.n are extras.
        idxs = idxs[idxs != tree.n]

        print sline.rstrip()
        for i in idxs:
            iline = info_lines[i]
            ikey = (iline[1], iline[5])
            if ikey in seen: continue
            seen[ikey] = True
            print "\t".join(iline)


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-d", "--dist", dest="dist", help="distance around each pair to look for missed pairs")
    parser.add_option("--diags", dest="diags", help="the dagchainter output file")
    parser.add_option("--all", dest="all", help="the dag blast hit file containing all hits of q to s")

    (options, _) = parser.parse_args()
    if not (options.dist and options.diags and options.all):
        sys.exit(parser.print_help())

    main(int(options.dist), options.diags, options.all)
