import collections
import sys
from __init__ import DagLine


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query   = args[0]
        self.subject = args[1]

        self.pctid  = float(args[2])
        self.hitlen = int(args[3])

        self.nmismatch = int(args[4])
        self.ngaps     = int(args[5])

        self.qstart = int(args[6])
        self.qstop  = int(args[7])
        self.sstart = int(args[8])
        self.sstop  = int(args[9])

        self.evalue = float(args[10])
        self.score  = float(args[11])

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % (self.query, self.subject, self.eval, self.score)
    def __str__(self):
        return "\t".join(map(str, (getattr(self, attr) for attr in BlastLine.__slots__)))

def adjust_evalue(afile, expected_count=8, evalue_cutoff=5, oclass=DagLine, 
                  out=sys.stdout
                 ):
    """
    adjust the evalues in a dagchainer/blast file by the number of times they occur.
    query/subjects that appear often will have the evalues raise (made less 
    significant).
    `afile`: path to the blast or dag file.
    `expected_count`: the number of blast hits expect per accn. a good estimate
                      is per file: [number of lines ] / [ number of unique accns ]
                      making the value higher results in less evalue adjustment.    
    `evalue_cutoff`: dont print lines with an adjusted evalue above this
    `oclass`: either DagLine or BlastLine
    `out`: a filehandle. defaults to stdout. if a string. it's opened.
    """
    if oclass is BlastLine:
        name1, name2 = ('query', 'subject')
    else:
        name1, name2 = ('a_accn', 'b_accn')

    if isinstance(out, basestring): out = open(out, 'w')

    expected_count = float(expected_count)
    counts = collections.defaultdict(int)
    for line in open(afile):
        b = oclass(line)
        counts[getattr(b, name1)] += 1
        counts[getattr(b, name2)] += 1
    for line in open(afile):
        b = oclass(line)
        count = counts[getattr(b, name1)] + counts[getattr(b, name2)]
        b.evalue = b.evalue ** (expected_count / count)

        if b.evalue < evalue_cutoff: print >>out, str(b)



if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser("""
    adjust the evalues in a dagchainer/blast file by the number of times they occur.
    query/subjects that appear often will have the evalues raise (made less 
    significant).
    """)

    p.add_option('-f', dest='afile', help="path to a blast or dag file")
    p.add_option('-t', dest='type', help="'blast' or 'dag'.\
       if not specified, the extention the extension of the file is used.",
            default=None)
    p.add_option('-c', dest='expected', type='float', help=
                "the number of blast hits expect per accn. a good estimate\n"
                "is per file: [number of lines ] / [ number of unique accns ]\n"
                "making the value higher results in less evalue adjustment.\n"
                "values between 4 and 15 seem to work well.", default=8)
    p.add_option('-o', dest='out', help="if specified, a file to write to."\
                                     " if not specified, stdout is used", 
                                    default=None)

    p.add_option('-e', dest='evalue', type='float', help="filter (dont print)"
        " lines with a value higher than this", default=5)

    opts, _ = p.parse_args()
    if not opts.afile: sys.exit(p.print_help())
    if not opts.type:
        if opts.afile.endswith(".blast"): opts.type = "blast"
        elif ".dag" in opts.afile: opts.type = "dag"
        else:
            sys.exit(p.print_help("specify the file type 'blast' or 'dag'"))
    assert opts.type in ("dag", "blast")
    if opts.type == "dag": oclass = DagLine
    else: oclass = BlastLine

    if opts.out is not None:
        out = opts.out
    else: 
        out = sys.stdout 
    adjust_evalue(opts.afile, opts.expected, opts.evalue, oclass, out)
