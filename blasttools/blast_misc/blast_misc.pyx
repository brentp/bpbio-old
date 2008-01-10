"""
take a tab-delimited blast hit and turn the relative coordinates into global,
absolute coordinates, the name format (columns 1, 2) must be in the
format: 
    chr||name||start||stop
    01||AT2G26540||123680||124220
"""

import operator 
import collections
import os
import numpy 
import cPickle


BLAST_NAMES = ['query', 'subject', 'percent_id', 'hit_length', 'n_mismatch', \
    'n_gaps', 'qstart', 'qstop', 'sstart', 'sstop', 'eval', 'score','hit_rank','dup']

BLAST_TYPES = ['S8', 'S8', 'f4', 'i4', 'i4', 'i4', 'i8', 'i8', 'i8', 'i8', 'f8', 'f8', 'i4', 'b1']


cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass
    FILE * fopen(char *, char *)
    int fscanf(FILE * i, char * fmt, ...) 
    cdef char * fgets(char *str, int s, FILE *f)
    cdef int fclose(FILE *)
    cdef Py_ssize_t strlen(char *)
    int EOF


cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass
    ctypedef struct PyListObject:
        pass
    ctypedef struct PyFloatObject:
        pass

    int PyList_Append(object l, object i) except -1
    object PyList_New(int)
    cdef int PyList_SET_ITEM(object, int, object) except -1
    cdef void Py_INCREF(void *)
    cdef object PyDict_GetItem(object d, object k)
    

cdef extern from "cblast_misc.c":
    int add_locs(char *, char *)

cdef char * blast_format = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%f*%s\n"

#cdef inline int set_list(object PyList, int idx, object obj) except -1:
#    Py_INCREF(<PyObject *>obj)
#    PyList_SET_ITEM(PyList, idx, obj)


def pyaddlocs(inname, outname=None):
    """ take a file name and add the global locs 
    inplace if outname==None, or to a new file (outname).
    """
    if outname is None:
        return add_locs(inname, NULL)
    return add_locs(inname, outname)

def pylocslist(*args):
    """ take a list of files to update inplace """
    if isinstance(args[0], list):
        args = args[0]
    for a in args:
        print add_locs(a, NULL)


def _blast_array(blastf, maxkeep):
    cdef int qlen = 8, slen = 8, hlen, nmiss, ngap, qstart, qstop, sstart, sstop, icol = 0
    cdef float pct = 0.0, evalue = 0.0, bit = 0.0
    cdef char qname[128], sname[128]
    cdef char *tmp
    cdef FILE *fh

    if maxkeep == None: maxkeep = 6

    recs = collections.defaultdict(list)
    fh = fopen(blastf, "r")

    while fscanf(fh, blast_format, qname, sname, &pct, &hlen, &nmiss, &ngap, &qstart, &qstop, &sstart, &sstop, &evalue, &bit ) != EOF:
        if strlen(qname) > qlen: qlen = strlen(qname)
        if strlen(sname) > slen: slen = strlen(sname)

        #key = qname + '_' + sname
        PyList_Append(recs[qname + '_' + sname], [qname, sname, pct, hlen, nmiss, ngap, qstart, qstop, sstart, sstop, evalue, bit, None, False])
    fclose(fh)

    arr = [] #PyList_New(0)
    # this is faster then using k,v in rec.items()
    for k in sorted(recs, key=operator.itemgetter(1)):
        
        v = PyDict_GetItem(recs, k)  # was: v = recs[k]
        icol = 0
        for li in sorted(v, key=operator.itemgetter(10))[:maxkeep]:
            # set the first to have hit_rank=0
            li[12] = icol
            icol += 1
            PyList_Append(arr, li)
    return arr, qlen, slen


def blast_array(fi, dopickle=1, best_hit=1, maxkeep=6):
    file_name = 0
    if not os.path.exists(fi): return []
    if not isinstance(fi, list) and type(fi) != numpy.recarray:
        file_name = fi + ".pickle"

        if dopickle and os.path.exists(file_name):
            if not best_hit:
                return numpy.load(file_name)
            else:
                arr = numpy.load(file_name)
                return arr[arr['hit_rank'] == 0]

        arr, qlen, slen = _blast_array(fi, maxkeep)
        BLAST_TYPES[0] = 'S' + str(qlen)
        BLAST_TYPES[1] = 'S' + str(slen)
    else:
        arr = fi
    if not len(arr):
        return []

    if type(arr) == numpy.recarray:
        ra = arr
    else:
        ra = numpy.rec.array(arr, names=BLAST_NAMES, formats=BLAST_TYPES)

    if dopickle and file_name: ra.dump(file_name)
    # hit_rank of 0 is the best hit
    if best_hit: ra = ra[ra['hit_rank'] == 0]
    return ra
 
def geneorder(fn, picklef=None, dups=set([]), field=None):
    """\
    take an array (`fn`) from blast_array and **returns the ordered genelist**
    as a hash where the values are the index/order.  this should be called with
    all self-self blast chromosomes to make sure to get all accns.
    **NOTE**: this assumes that the genes can be sorted lexically as in
    arabidopsis: ['AT1G01010','AT1G01020' ...]
    if they cannot be sorted lexically, you must provide a new function which
    returns a dict of name=>order pairs.
    """
    arr = blast_array(fn, best_hit=1)
    if field:
        genelist = numpy.sort(numpy.unique( arr.field(field)))
    else:
        genelist = numpy.sort(numpy.unique(numpy.vstack(( arr.field("query"), arr.field("subject") ))   ))
    gorder = dict([(g,i) for g, i in zip(genelist, xrange(genelist.shape[0])) \
                      if not g in dups])

    if picklef is not None:
        cPickle.dump(gorder, open(picklef, 'wb'), -1)
    return gorder


def order_from_blast_list(*args):
    """ args is a list of blast filenames
    not yet implemented """
    cdef int qlen = 8, slen = 8, hlen, nmiss, ngap, qstart, qstop, sstart, sstop, icol = 0
    cdef float pct = 0.0, evalue = 0.0, bit = 0.0
    cdef char qname[128], sname[128]
    cdef FILE *fh
    for f in args:
        fh = fopen(f, "r")
        while fscanf(fh, blast_format, qname, sname, &pct, &hlen, &nmiss, &ngap, &qstart, &qstop, &sstart, &sstop, &evalue, &bit ) != EOF:
           pass 
                    
        fclose(fh)
