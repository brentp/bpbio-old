"""
take a genomic `fasta` and the `blast` output which resulted from the self-self
blast of that fasta. and generate a new fasta in which all basepairs occuring
more then `cutoff` times are masked. e.g.:
    $ python mask_genome.py -b rice_rice.blast -o rice -f rice.fasta -c 50
will create a new fasta file rice.masked.50.fasta
""" 
import numpy as np
from pyfasta import Fasta
import tables
import numexpr
import sys, os

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

H5 = 'copy_count.h5'
def get_node(org, mode):
    # get the parent group node in the h5 file.
    if mode == 'w':
        h5 = tables.openFile(H5, mode='w')
        if org in h5.root:
            raise Exception('that organism already exists in %s' % H5)
        return h5, h5.createGroup(h5.root, org, org)
    else:
        h5 = tables.openFile(H5, mode='r')
        return h5, getattr(h5.root, org)

def count_freq(blast_file, fasta, org):
    """one large blast file """
    h5, node = get_node(org, 'w')
    f = Fasta(fasta)

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

    h5.close()

def mask(fasta_file, org, cutoff, mask_value='X'):
    h5, node = get_node(org, 'r')

    outfile = fasta_file[:fasta_file.rfind(".")] + (".masked.%i" % cutoff) \
                         + fasta_file[fasta_file.rfind("."):]

    print "creating file:", outfile
    out = open(outfile ,'w')

    fasta = Fasta(fasta_file)

    soft_mask = mask_value.lower() == 'soft'
    for seqid in sorted(fasta.iterkeys()): 

        if soft_mask:
            seq = str(fasta[seqid])
            # mask is the lowercase sequence.
            mask_value = np.array(seq.lower(), dtype='c')
            seq = np.array(seq.upper(), dtype='c')
        else:
            fasta[seqid].tostring = False
            seq = fasta[seqid][:] # a


        if not 'c' + seqid in node:
            print >>sys.stderr, seqid, 'not found in masked, writing unchanged'
            out.write('>' + seqid + '\n')
            out.write(seq.tostring() + '\n')
            continue
        
        hit_counts = getattr(node, 'c' + seqid)[:]
        masked_seq = np.where(numexpr.evaluate("hit_counts > %i" % cutoff)
                              , mask_value, seq).tostring() 

        print >>sys.stderr, seqid, len(seq)
        assert len(seq) == len(masked_seq)
        out.write('>' + seqid + '\n')
        out.write(masked_seq + '\n')

    out.close()
    h5.close()

if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("-b", dest="blast", help="path to self-self genomic blast")
    p.add_option("-f", dest="fasta", help="path to the fasta file (which was\n"
                 "used as query and subject in the blast")
    p.add_option("-o", dest="org", help="name of the organism. e.g. 'rice'")
    p.add_option("-c", dest="cutoff", 
                 help="cutoff value, bp locations appearing in the blast more\n"
                 "than this many times are masked", type='int', default=50)
    p.add_option("-m", dest="mask", help=\
         "mask sequence with this letter (usually 'X' or 'N'). if == 'SOFT',"
         "then the sequence subjected to soft-masking where all repetitive"
         "values are lower-cased. and all other sequence is upper-cased"
         "regardless of its case in the original fasta", default='X')

    options, _ = p.parse_args()
    if not (options.blast and options.fasta and options.org):
        sys.exit(p.print_help())
    ospe = os.path.exists
    if not (ospe(options.blast) and ospe(options.fasta)):
            print "make sure blast:%s and fasta:%s exist:" \
                                % (options.blast, options.fasta)
            sys.exit()
    assert len(options.mask) == 1 or options.mask.lower() == 'soft'
    count_freq(options.blast, options.fasta, options.org)
    print "DONE COUNTING"
    mask(options.fasta, options.org, options.cutoff, options.mask)
