import numpy as np
from pyfasta import Fasta
import tables
import numexpr

def cache_clear(cache, node, qchr, schr):
    """ keep some of the arrays in memory
    if there are too many, put one back into
    the h5 file
    """
    if len(cache) < 5: return
    #print "updating cache: %s, %s" % (qchr, schr)
    rmkey = [k for k in cache if not k in (qchr, schr)][0]
    getattr(node, 'c' + rmkey)[:] = cache[rmkey]
    del cache[rmkey]

def update_cache(achr, node, clen, h5, cache):
    """
    put a new array in memory, and make sure it's
    in the h5 as well
    """
    if not 'c' + achr in node:
        h5.createArray(node, 'c' + achr, np.zeros((clen,), dtype=np.uint32))

    cache[achr] = getattr(node, 'c' + achr)[:]


def count_freq(blast_file, fasta):
    """one large blast file """
    h5 = tables.openFile(blast_file + '.h5', mode='w')
    h5.createGroup(h5.root, 'count', 'count')
    f = Fasta(fasta)

    node = h5.root._v_children['count']

    cache = {}
    for sline in open(blast_file):
        line = sline.split("\t")
        qchr, schr = line[:2]

        qstart, qstop, sstart, sstop = map(int, line[6:10])
        if sstart > sstop: sstart, sstop = sstop, sstart

        if not qchr in cache:
            update_cache(qchr, node, len(f[qchr]), h5, cache)
            cache_clear(cache, node, qchr, schr)

        if not schr in cache:
            update_cache(schr, node, len(f[schr]), h5, cache)
            cache_clear(cache, node, qchr, schr)

        cache[qchr][qstart - 1: qstop - 1] += 1
        cache[schr][sstart - 1: sstop - 1] += 1

    for achr in cache:
        getattr(node, 'c' + achr)[:] = cache[achr]

    return node

def mask(node, fasta_file, org, cutoff=50):


    outfile = fasta_file[:fasta_file.rfind(".")] + ".masked" \
                         + fasta_file[fasta_file.rfind("."):]

    print "creating file:", outfile
    out = open(outfile ,'w')

    fasta = Fasta(fasta_file)


    for seqid in sorted(fasta.iterkeys()): 
        fasta[seqid].tostring = False
        seq = fasta[seqid][:] # a

        if not 'c' + seqid in node:
            print >>sys.stderr, seqid, 'not found in masked, writing unchanged'
            out.write('>' + seqid + '\n')
            out.write(seq.tostring() + '\n')
            continue
        
        hit_counts = getattr(node, 'c' + seqid)[:]
        masked_seq = np.where(numexpr.evaluate("hit_counts > %i" % cutoff), "X", seq).tostring() # + add
        print seqid, len(seq)
        assert len(seq) == len(masked_seq)
        out.write('>' + seqid + '\n')
        out.write(masked_seq + '\n')

    out.close()

if __name__ == "__main__":
    import sys
    # blast, fasta, org
    blast, fasta, org = sys.argv[1], sys.argv[2], sys.argv[3]
    node = count_freq(blast, fasta)
    print "DONE COUNTING"
    mask(node, fasta, org)
