r"""

    >>> b = BlastLine('chr6\t10\t91.81\t1758\t144\t0\t115531278\t115533035\t9465094\t9466851\t0.0\t2343.0')
    >>> b
    BlastLine(chr6[115531278:115533035]-10[9465094:9466851], ptcid=91.810, eval=0.000)

"""
from feature import *
#
class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query  =args[0]
        self.subject  = args[1]
        self.pctid =float(args[2])
        self.hitlen =int(args[3])
        self.nmismatch =int(args[4])
        self.ngaps =int(args[5])
        self.qstart =int(args[6])
        self.qstop =int(args[7])
        self.sstart =int(args[8])
        self.sstop =int(args[9])
        self.evalue =float(args[10])
        self.score =float(args[11])

    def __repr__(self):
        return ("BlastLine(%s[%i:%i]-%s[%i:%i], " +\
                   "ptcid=%.3f, eval=%.3f)") \
            % (self.query, self.qstart, self.qstop, self.subject, self.sstart,\
               self.sstop, self.pctid,
               self.evalue)

try:
    from biostuff import BlastLine
except:
    pass

