"""
take any number of non-sorted gff3 files and spit out a sorted
version that genometools will accept.
my need to add more top-level features to sortable_type()
"""
import sys
from biostuff import GFFLine

files = sys.argv[1:]

lines = {}

def sort_feats(a, b):
    if a.type in ('gene', 'pseudogene', 'MIR'): return -2
    if b.type in ('gene', 'pseudogene', 'MIR'): return 2
    if a.attribs.get('ID') and not b.attribs.get('ID'): return - 1
    if b.attribs.get('ID') and not a.attribs.get('ID'): return 1
    if a.type == 'mRNA' and b.type in ('mRNA', 'exon', 'CDS'): return -1
    if b.type == 'mRNA' and a.type in ('mRNA', 'exon', 'CDS'): return 1
    return cmp(a.start, b.start) or cmp(a.end, b.end)

def find_exts(flist):
    fmin = min(f.start for f in flist)
    fmax = max(f.end for f in flist)
    return fmin, fmax

for fi in files:
    if fi == "-": fi = sys.stdin
    else: fi = open(fi)
    for line in fi:
        if line[0] == '#': continue
        o = GFFLine(line)
        fid = o.attribs.get('Parent', o.attribs.get('rname', o.attribs.get('ID')))
        # hack...
        while 'mRNA' in fid and fid[-1].isdigit(): # something.mRNA.1
            fid = fid[:-2]

        fid = fid.rstrip('.mRNA')

        if fid is None:
            raise Exception("No ID %s" % str(o))
        if not o.seqid in lines:
            lines[o.seqid] = {}
        if not fid in lines[o.seqid]:
            lines[o.seqid][fid] = []
        lines[o.seqid][fid].append(o)

def fcmp(a, b):
    return cmp(min(f.start for f in a[1]), min(f.start for f in b[1]))

print '##gff-version 3'

seen = {}
for seqid, iddicts in sorted(lines.iteritems()):
    for nameid, features in sorted(iddicts.items(), cmp=fcmp):

        fmin, fmax = find_exts(features)
        for f in features:
            if f.type == 'gene':
                if fmin < f.start: f.start = fmin
                if fmax > f.end:   f.end   = fmax

        for f in sorted(features, cmp=sort_feats):
            # 9/30/09 added ID because sometimes have same gene with diff id?
            key = (f.seqid, f.start, f.end, f.type, f.attrs.get("ID"))
            if key in seen: 
                print >>sys.stderr, key
                continue
            seen[key] = None
            print f.to_line()
