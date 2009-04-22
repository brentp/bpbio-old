from cblastline import BlastLine
from rtree import Rtree
from gff_reader import GFFDict


def blast_to_tree(blast_file, qgff=None, sgff=None):
    """
    create a series of rtree's from a blast file
    the gff files are used to adjust the blast positions
    from local to global.
    if the blast position are chromosome based, the gff files
    should not be specified.
    """

    r = {}
    counts = {}
    if qgff:
        qgff = GFFDict(qgff)
    if sgff:
        sgff = GFFDict(sgff)

    for line in open(blast_file):
        b = BlastLine(line)

        if not b.query in r:
            r[b.query] = {}
            counts[b.query] = {}
        if not b.subject in r[b.query]:
            r[b.query][b.subject] = {'blasts':[], 'tree': Rtree()}
            counts[b.query][b.subject] = 0
        d = r[b.query][b.subject]
        d['blasts'].append(b)

        smin, smax = b.sstart, b.sstop
        if smin > smax: smax, smin = smin, smax

        qmin, qmax = b.qstart, b.qstop

        if sgff is not None:
            s = sgff[b.subject]
            smin += s.start
            smax += s.start
        if qgff is not None:
            q = qgff[b.query]
            qmin += q.start
            qmax += q.start

        #r[b.query][b.subject]['tree'].add(counts[b.query][b.subject], (b.qstart, smin, b.qstop, smax))
        d['tree'].add(counts[b.query][b.subject], (qmin, smin, qmax, smax))
        counts[b.query][b.subject] += 1
    return r


