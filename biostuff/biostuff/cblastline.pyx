cdef extern from *:
    ctypedef char* const_char_star "const char*"

cdef extern from "stdio.h":
    int sscanf(char* astr, const_char_star format, ...)

cdef extern from "Python.h":
    char *PyString_AsString(object)

cdef const_char_star blast_format = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%f"



cdef class BlastLine:
    r"""
    given a string of tab-delimited (-m 8) blast output, parse it and create
    an object with the usual attrs:

    >>> from biostuff import BlastLine 
    >>> b = BlastLine("Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0")
    >>> b.query
    'Os09g11510'

    >>> attrs = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
    ...  'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')
    >>> [getattr(b, attr) for attr in attrs]# doctest: +ELLIPSIS
    ['Os09g11510', 'Os08g13650', 92.3..., 39, 3, 0, 2273, 2311, 3237, 3199, 0.001..., 54.0]
        
    """
    cdef public int hitlen, nmismatch, ngaps, qstart, qstop, sstart, sstop
    cdef public float pctid, score, evalue
    cdef char _cquery[48], _csubject[48]
    cdef object _pysubject, _pyquery

    property query:
        def __get__(self):
            if self._pyquery is None:
                return self._cquery
            return self._pyquery
        def __set__(self, val):
            self._pyquery = val
    property subject:
        def __get__(self):
            if self._pysubject is None:
                return self._csubject
            return self._pysubject
        def __set__(self, val):
            self._pysubject = val

    def __init__(self, char *sline):
        sscanf(sline, blast_format, self._cquery, self._csubject,
                &self.pctid, &self.hitlen, &self.nmismatch, &self.ngaps,
                &self.qstart, &self.qstop,
                &self.sstart, &self.sstop,
                &self.evalue, &self.score)

    def __repr__(self):
        return "BlastLine('%s'-'%s', ptcid=%.3f, eval=%.3f, score=%.1f)" \
            % (self.query, self.subject, self.pctid, self.evalue, self.score)


    def to_blast_line(self, as_str=True):
        attrs = ['query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score']
        if as_str:
            return "\t".join(map(str, [getattr(self, attr) for attr in attrs]))
        else:
            return [getattr(self, attr) for attr in attrs]

