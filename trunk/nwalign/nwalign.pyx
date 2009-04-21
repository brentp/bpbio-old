import numpy as np
cimport numpy as np


cimport cython
import sys


cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void free(void *ptr)

cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *, size_t)

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT

cdef size_t UP = 0, LEFT = 1, DIAG = 2, NONE = 3

MATRIX = {}

cdef inline int imax2(int a, int b):
    if a >= b: return a
    return b


def read_matrix(path):
    if path in MATRIX: return MATRIX[path]
    m = {}
    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [x for x in line.split(' ') if x]
        for h in headers: m[h] = {}

    line = fh.readline()
    while line:
        h1 = line[0]
        line = [int(x) for x in line[1:-1].split(' ') if x]
        values = zip(headers, line)
        m[h1] = dict(values)
        line = fh.readline()
    return m

@cython.boundscheck(False)
cdef global_align(object _seqj, object _seqi, int gap, int match, int mismatch, object matrix):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')
    """
    cdef char* seqj = _seqj
    cdef char* seqi = _seqi

    cdef int max_j = strlen(seqj)
    cdef int max_i = strlen(seqi)
    cdef size_t i, j, seqlen, align_counter = 0, p
    cdef int diag_score, up_score, left_score, tscore

    cdef char *align_j, *align_i
    cdef char ci, cj


    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)


    if matrix is not None:
        matrix = read_matrix(matrix)
  
    #score   = np.zeros((max_i + 1, max_j + 1), dtype='f')
    #pointer = np.zeros((max_i + 1, max_j + 1), dtype='i')

    pointer[<size_t>0, <size_t>0] = NONE
    score[<size_t>0, <size_t>0] = 0

    
    pointer[<size_t>0, <size_t>1:] = LEFT
    pointer[<size_t>1:, <size_t>0] = UP

    score[<size_t>0, <size_t>1:] = gap * np.arange(max_j, dtype=np.int)
    score[<size_t>1:, <size_t>0] = gap * np.arange(max_i, dtype=np.int)
    
    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]

            if matrix is None:
                diag_score = score[i - 1, j - 1] + (cj == ci and match or mismatch)
            else:
                try:
                    py_ci = PyString_FromStringAndSize(&ci, <size_t>1)
                    py_cj = PyString_FromStringAndSize(&cj, <size_t>1)
                    tscore = matrix[py_cj][py_ci]
                    diag_score = score[i - 1, j - 1] + tscore
                except:
                    print (cj + c'0'), ci, "not in matrix"
                    raise

            up_score   = score[i - 1, j] + gap
            left_score = score[i, j - 1] + gap
            
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                else:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT

            else:
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
    seqlen = imax2(i, j)
    align_j = <char *>malloc(seqlen * sizeof(char))
    align_i = <char *>malloc(seqlen * sizeof(char))
        
    while True:
        p = pointer[i, j]
        if p == NONE: break
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = c"-"
        elif p == UP:
            i -= 1
            align_j[align_counter] = c"-"
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!')
        align_counter += 1

    try:
        return str(align_j)[:align_counter][::-1], str(align_i)[:align_counter][::-1]
    finally:
        free(align_j)
        free(align_i)
            
def main():
    import optparse
    parser = optparse.OptionParser(usage="""
    %prog [options] seq1 seq2 
    """)
    parser.add_option("--gap", dest="gap", help="gap penalty (must be integer < 0)", type="int", default=-1)
    parser.add_option("--match", dest="match", help="match score (must be integer > 0)", type="int", default=1)
    parser.add_option("--mismatch", dest="mismatch", help="gap penalty (must be integer < 0)", type="int", default=-1)
    parser.add_option("--matrix", dest="matrix", help="scoring matrix in ncbi/data/ format,\
                                      if not specificied, match/mismatch are used", default=None)

    try:
        options, args = parser.parse_args()
    except:
        sys.exit(parser.print_help())

    if len(args) != 2:
        sys.exit(parser.print_help())
    print "\n".join(global_align(args[0], args[1], options.gap, options.match, options.mismatch, options.matrix))
        
