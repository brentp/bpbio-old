cdef extern from *:
    ctypedef char* const_char_star "const char*"

cimport stdlib

cdef const_char_star dag_format_line = "%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%lf"

cdef extern from "stdio.h":
    int sscanf(char* astr, const_char_star format, ...)

DEF MAXSIZE=64

cdef class DagLine:
    r"""
    given a string of tab-delimited dag output, parse it and create
    an object with the usual attrs.
    """
    cdef public int a_start, a_end, b_start, b_end
    cdef public double evalue
    cdef public char *a_accn, *a_seqid, *b_accn, *b_seqid

    def __init__(self, char *sline=NULL):
        cdef char a_accn[MAXSIZE]
        cdef char b_accn[MAXSIZE]
        cdef char a_seqid[MAXSIZE]
        cdef char b_seqid[MAXSIZE]

        if sline == NULL: return
        sscanf(sline, dag_format_line, 
               a_seqid, a_accn, &self.a_start, &self.a_end,
               b_seqid, b_accn, &self.b_start, &self.b_end,
               &self.evalue)
        self.a_accn = <char *>(stdlib.malloc(sizeof(char *) * stdlib.strlen(a_accn)))
        self.b_accn = <char *>(stdlib.malloc(sizeof(char *) * stdlib.strlen(b_accn)))
        self.a_seqid = <char *>(stdlib.malloc(sizeof(char *) * stdlib.strlen(a_seqid)))
        self.b_seqid = <char *>(stdlib.malloc(sizeof(char *) * stdlib.strlen(b_seqid)))
        stdlib.strcpy(self.a_accn, a_accn)
        stdlib.strcpy(self.b_accn, b_accn)
        stdlib.strcpy(self.a_seqid, a_seqid)
        stdlib.strcpy(self.b_seqid, b_seqid)
        if self.evalue < 1e-250: self.evalue = 1e-250


    def __dealloc__(self):
        stdlib.free(<void *>self.a_accn)
        stdlib.free(<void *>self.b_accn)
        stdlib.free(<void *>self.a_seqid)
        stdlib.free(<void *>self.b_seqid)

    def __repr__(self):
        return ("DagLine('%s', '%s')" % (self.a_accn, self.b_accn))

    @classmethod
    def from_dict(cls, dict d):
        return _factory(d)

cdef DagLine _factory(dict d):
    cdef DagLine instance = DagLine.__new__(DagLine)
    instance.a_seqid = <char *>stdlib.malloc(sizeof(char *) * len(d['a_seqid']))
    stdlib.strcpy(instance.a_seqid, d['a_seqid'])
    instance.b_seqid = <char *>stdlib.malloc(sizeof(char *) * len(d['b_seqid']))
    stdlib.strcpy(instance.b_seqid, d['b_seqid'])

    instance.a_accn = <char *>stdlib.malloc(sizeof(char *) * len(d['a_accn']))
    stdlib.strcpy(instance.a_accn, d['a_accn'])
    instance.b_accn = <char *>stdlib.malloc(sizeof(char *) * len(d['b_accn']))
    stdlib.strcpy(instance.b_accn, d['b_accn'])

    instance.a_start = d['a_start']
    instance.b_start = d['b_start']
    instance.a_end = d['a_end']
    instance.b_end = d['b_end']
    instance.evalue = d['evalue']
    return instance

