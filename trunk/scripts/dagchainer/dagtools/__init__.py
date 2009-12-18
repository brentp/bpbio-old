try:
    from cdagline import DagLine
except ImportError:
    # python version in case the cython one isnt compile.d
    class DagLine(object):
        __slots__ = ('a_seqid', 'a_accn', 'a_start', 'a_end', 
                     'b_seqid', 'b_accn', 'b_start', 'b_end', 'evalue')

        def __init__(self, line):

            tline = line[:-1].split("\t")
            self.a_seqid = tline[0]
            self.b_seqid = tline[4]
            self.a_accn = tline[1]
            self.b_accn = tline[5]

            self.a_start, self.a_end = map(int, tline[2:4])
            self.b_start, self.b_end = map(int, tline[6:8])
            # so if it doesnt have an evalue or it's empty,
            # just assume the best score.
            if len(tline) < 9 or not tline[8]: self.evalue = 1e-250
            else: self.evalue = max(float(tline[8]), 1e-250)
        def __repr__(self):
            return "DagLine('%s', '%s')" % (self.a_accn, self.b_accn)

        @classmethod
        def from_dict(cls, d):
            for k in cls.__slots__:
                setattr(cls, k, d[k])

