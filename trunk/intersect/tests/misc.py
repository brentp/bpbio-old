import collections

def read_gff(gff_file, chrn=None, name_fn=lambda l: l[8]):
    """\
    a function that takes a gff file, parses it, and returns a Genedex-able object of the
    requested chromosome.
    gff_file: ...

    chrn: and identifier for chromosome that matches the first column in gff_file, if None
          return all data.

    name_fn: a function that accepts the line array (after being split on tabs) and returns 
             either the name to use for that feature, or None, if
             None, the line is skipped. the default just uses the 8th
             column as the name. 
                >>> name_fn = lambda l: l[8]
             To use a regular expression to extract the name from the last column.
             e.g. if the last column is "Parent=AT1G01010;Note=protein_coding_gene;" a
             reasonable name_fn would be:
                >>> def name_fn(line):
                ...     match = re.search("=([^.;]+)", l[8], re.I)
                ...     if match is None: return None
                ...     else: return match.groups(0)[0]
            which would return "AT1G01010" 
    """
    gff_by_chr = collections.defaultdict(list)
    by_name = {}
    types = set([])
    for line in open(gff_file):
        line = line[:-1].split('\t')
        ftype = line[2].upper()
        types.update([ftype])
        feat = dict(
                start = int(line[3])
                , stop  = int(line[4])
                , type  = ftype
                , strand= line[6] in ('-', '-1') and -1 or 1
                , chr   = line[0]
                )
        feat[ftype] = [(feat['start'], feat['stop'])]
        n = name_fn(line)
        if n is None: continue
        feat['name'] =  n
        # so in a gff file, each exon/cds, etc is stored in a separate
        # line. here, it merges these by using the same feature and
        # appending the locs to the appropriate feature type.
        if feat['name'] in by_name:
            existing = by_name[feat['name']]
            if feat['type'] in existing:
                existing[feat['type']].extend(feat[ftype])
            else:
                existing[feat['type']] = feat[ftype]
        else:
            by_name[feat['name']] = feat

    # find the overall start, stop
    for feature in by_name.itervalues():
        for ftype in types:
            if not ftype in feature: continue
            locs = feature[ftype]
            start = min(x[0] for x in locs)
            stop  = max(x[1] for x in locs)
            if start < feature['start']: feature['start'] = start
            if stop  > feature['stop'] : feature['stop']  = stop

        gff_by_chr[feature['chr']].append(feature)

    if not chrn: return gff_by_chr
    return gff_by_chr[chrn]

