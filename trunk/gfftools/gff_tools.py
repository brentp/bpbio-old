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

