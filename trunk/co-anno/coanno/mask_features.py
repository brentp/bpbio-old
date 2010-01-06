from pyfasta import Fasta
import os
import sys
import numpy as np
import gtpym

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(os.path.basename(__file__))


def is_up_to_date(m, g, f):
    """
    >>> is_up_to_date('a.1.t', 'a.2.t', 'a.3.t')
    False
    >>> is_up_to_date(__file__, __file__, __file__)
    False
    """
    ospe = os.path.exists
    if not ospe(m): return False
    if not (ospe(g) and ospe(f)): return False
    mtime = os.stat(m).st_mtime
    return mtime > os.stat(g).st_mtime \
            and mtime > os.stat(f).st_mtime

def get_fastas(fasta_file, mask_types):
    """
    >>> get_fastas('a.fasta', ('CDS', 'gene'))
    ('a.genomic.masked.fasta', 'a.features.fasta')

    >>> get_fastas('a.fasta', [])
    ('a.genomic.masked.genes.fasta', False)
    """

    if mask_types:
        genomic_masked = fasta_file.replace('.fa', '.genomic.masked.fa')
        features_fasta = fasta_file.replace('.fa', '.features.fa')
    else:
        genomic_masked = fasta_file.replace('.fa', '.genomic.masked.genes.fa')
        features_fasta = False #fasta_file.replace('.fa', '.genes.fa')

    return genomic_masked, features_fasta

def write_feature(fh, feat, fasta, seen_ids, mask_types):
    r"""
    >>> from cStringIO import StringIO
    >>> fh = StringIO()
    >>> feat =  gtpym.FeatureNode('chr2', 'CDS', 2, 4, '+')
    >>> feat.add_attribute('ID', 'some_id')
    >>> fasta = Fasta('tests/data/t.fasta')
    >>> write_feature(fh, feat, fasta, {}, 'CDS')
    >>> fh.getvalue()
    '>some_id\nCTC\n'

    """
    fname = feat.get_attribute("ID") or feat.get_attribute('rname')
    assert not fname in seen_ids, (fname, "used more than 1X", feat)
    assert fname is not None, (feat, feat.type, feat.attribs )
    seen_ids[feat] = None
    fh.write('>' + fname + '\n')
    fdict = {'chr': feat.seqid, 'start': feat.start, 'stop': feat.end, 'strand': feat.strand }
    fh.write(fasta.sequence(fdict , exon_keys=mask_types) + "\n")


def check_exists(f, raw_input=raw_input):
    """
    >>> raw_input = lambda a: 'yes'
    >>> check_exists('a.ttttttttt')
    >>> check_exists(__file__, raw_input=raw_input)
    """
    if f and os.path.exists(f):
        r = raw_input("%s exists, do you wish to overwrite? (y/(n)) " % (f,))
        if r.lower()[0] != 'y':
            raise Exception("wont overwrite existing file without permission")

class NotMaskedException(Exception): pass

def mask_seq_with_types(sequence, mask_types, feat, N):
    """
    >>> s = np.array('AAAAAAAAAAAAAAAAAAAAAA', dtype='c')
    >>> N = np.array('N', dtype='|S1')
    >>> mask_types = ['CDS']
    >>> feat = gtpym.FeatureNode('chr2', 'CDS', 4, 8, '+')
    >>> mask_seq_with_types(s, mask_types, feat, N)
    4L

    >>> s.tostring()
    'AAANNNNNAAAAAAAAAAAAAA'

    >>> s = np.array('AAAAAAAAAAAAAAAAAAAAAA', dtype='c')
    >>> mask_seq_with_types(s, None, feat, N)
    4L

    >>> s.tostring()
    'AAANNNNNAAAAAAAAAAAAAA'

    """
    feats_by_type = gtpym.FeatureIndexMemory.by_types(feat)

    n_masked = 0
    if mask_types:
        for ftype in mask_types:
            if not ftype in feats_by_type: continue
            locs = feats_by_type[ftype]

            for loc in locs:
                sequence[loc.start - 1: loc.end] = N
                n_masked += (loc.end - loc.start)
            return n_masked

        else: # not masked
            raise NotMaskedException
    else:
        # mask the whole thing.
        sequence[feat.start - 1: feat.end] = N
        return (feat.end - feat.start)


def main(gff_file, fasta_file, mask_types=('CDS', 'mRNA', 'MIR', 'gene'),
                                            inverse=False):
    genomic_masked, features_fasta = get_fastas(fasta_file, mask_types)

    if is_up_to_date(genomic_masked, gff_file, fasta_file) \
            and is_up_to_date(features_fasta, gff_file, fasta_file):
        log.debug("%s is up-to-date." % (genomic_masked, ))
        return False

    gff = gtpym.FeatureIndexMemory(gff_file)
    fasta = Fasta(fasta_file)


    N = np.array('@' if inverse else 'N', dtype='c')
    if inverse:
        genomic_masked = genomic_masked.replace('genomic.masked', 'genomic.nonfeat.masked')

    seen_ids = {}
    for f in (genomic_masked, features_fasta):
        check_exists(f)

    try:
        genomic_masked  = open(genomic_masked, "w")
        features_fasta = f and open(features_fasta, "w")

        fkeys = sorted(fasta.keys())
        if len(fkeys) >= 100:
            log.debug("beginning masking for %i sequences" % (len(fkeys,)))

        for achr in fkeys:
            features = gff.get_features_for_seqid(achr)
            sequence = np.array(fasta[achr])

            tot_masked = 0
            for feat in features:
                try:
                    tot_masked += mask_seq_with_types(sequence, mask_types, feat, N)
                except NotMaskedException:
                    print >>sys.stderr, feat
                    cleanup(genomic_masked, features_fasta)
                    raise NotMaskedException("the sequence %s was not masked" \
                                             % (feat.get_attribute("ID"), ))


                if features_fasta:
                    write_feature(features_fasta, feat, fasta, seen_ids, mask_types)

            if inverse:
                nseq = np.array(fasta[achr])
                nseq[sequence != '@'] = np.array('N', dtype='c')
                sequence = nseq.copy()
                assert sequence[sequence != 'N'].shape[0] > 0, (achr, sequence[:100])

            genomic_masked.write('>' + achr + "\n")
            genomic_masked.write(sequence.tostring() + "\n")
            if len(fkeys) < 100:
                log.debug("%i basepairs of %i (%.1f%%) masked for %s" % \
                        (tot_masked, sequence.shape[0],
                            100. * float(tot_masked)/sequence.shape[0], achr))
    except:
        cleanup(genomic_masked, features_fasta)
        raise


def cleanup(*files):
    """
    >>> fh = open('ttt.test', 'w')
    >>> fh.close()
    >>> cleanup(fh)
    """
    for f in files:
        if hasattr(f, 'name'): f = f.name
        if os.path.exists(f): os.unlink(f)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print >>sys.stderr, """\
   usage: %s path/to.gff path/to.fasta [mask feature_types] --inverse

        feature types can be a list of names like: CDS mRNA gene
        in which case, the first of those that is found will be used.
        to mask the sequence in the fasta file. 
        if the list is empty, the entire range of each feature is masked.
        if the last argument is '--inverse', everythign _but_ the feature types
        is masked.
                """ % sys.argv[0]
        sys.exit()

    gff = sys.argv[1].rstrip()
    fasta = sys.argv[2].rstrip()

    inverse = sys.argv[-1] == '--inverse'
    if len(sys.argv) > 3:
        mask_types = [x.strip() for x in sys.argv[3:]]
        main(gff, fasta, mask_types, inverse)
    else:
        main(gff, fasta, None, inverse)
