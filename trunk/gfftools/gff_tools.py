import collections
import operator

def get_sequence(feat, fasta, types=('CDS', 'mRNA', 'gene')):
    """
    get the genomic sequence associated with a feature.

    `feat` is a genometools (gt) feat.
    `fasta` is a pyfasta Fasta object
    `types` is an ordered list of feature types. as soon as the first one is 
            found in the feature, subfeatures of that type will be sorted by 
            start and used to create a single string sequence.
    """

    by_type = collections.defaultdict(list)
    for f in feat(method='depth_first'):
        by_type[f.type].append(f)
    by_type = dict(by_type)
    if by_type == {}:
        by_type = {'gene': [feat]}
    fa = fasta[feat.seqid]
    seq = ""
    for type in types:
        if type in by_type:
            fs = sorted(by_type[type], key=operator.attrgetter('start'))
            for f in fs:
                seq += fa[f.start - 1: f.end]
            return seq.upper()
    raise Exception("no sequence for %s" % feat.attribs["ID"])

from pyfasta import complement
import collections

def locs_for_feature(parent, prefs=('CDS', 'mRNA', 'gene', 'MIR')):
    try:
        types = collections.defaultdict(list)
        for f in parent:
            types[f.type].append(f)
        types = dict(types)

        for p in prefs:
            if p in types:
                return types[p], p
        else:
            return [parent], getattr(parent, 'type', 'type')
    except:
        return [parent], getattr(parent, 'type', 'type')


def sequence_for_feature(f, fa, locs=None):
    seq = ""

    fas = fa[f.seqid]

    if locs is None:
        locs, type = locs_for_feature(f)
    for loc in locs:
        seq += fas[loc.start -1: loc.end]
    if getattr(f, 'strand', '+') in ('-', -1, '-1'):
        return complement(seq)[::-1]
    return seq

def sequence_with_locs(seq, locs, as_idx=True, strand='+'):
    """
    return the global positions zipped with:
        the sequence when `as_idx` is False
        the local positoins when `as_idx` is True
    useful to see where the gaps/introns are"""
    
    locslist = []
    for l in locs:
        locslist += range(l.start - 1, l.end)
    assert len(locslist) == len(seq), (len(locslist), len(seq))
    if strand == '-':
        locslist = locslist[::-1]
    else: assert strand == '+'

    if as_idx:
        return list(enumerate(locslist))
    else:
        return zip(seq, locslist)

