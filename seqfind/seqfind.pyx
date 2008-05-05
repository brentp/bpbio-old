"""
This is a fast, cython implementation of this:
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/572156

to allow fast searching for similar sequences. 

provides:
    function edit_distance: which takes 2 strings, and a
                maximum distance (int) and returns the damerau-levenshtein distance
                between those strings:

    class BKtree: takes a list of ["words", "and", "creates", "a", "search", "tree"]
                  so that it can be searched for any words within dist
                  of a given word:
            >>> tree.find("asd", 1)
            ["and"]
            >>> tree.find("asd", 2)
            ["and", "a"]
"""

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b) 

cdef extern from "Python.h":
    void PyList_Append(object, object)
    object PyTuple_GET_ITEM(object, int)
    void Py_INCREF(object)

cdef inline size_t imin(int a, int b, int c):
    if a < b:
        if c < a:
            return c
        return a
    if c < b:
        return c
    return b


cpdef int edit_distance(char *a, char *b, int limit):
    """
    the damerau-levenshtein distance:
    http://en.wikipedia.org/wiki/Damerau-Levenshtein_distance
    which differs from levenshtein in having only 1 point for a
    transposition:
        "zzabzz" => "zzbazz" levenshtein == 2
        "zzabzz" => "zzbazz" damerau-lev == 1

    takes 2 strings and a limit, if the strings have a length
    difference > limit

    >>> edit_distance("asdf", "acdf", 100)
    1

    >>> edit_distance("acdf", "acdf", 100)
    0

    >>> edit_distance("ac", "ca", 100)
    1
    """



    cdef int alen = strlen(a), blen = strlen(b)
    cdef char *ctmp
    cdef size_t i, j, retval
    cdef size_t achr, bchr #, cost
     
    if strcmp(a, b) == 0:
        return 0

    if alen > blen:
        ctmp = a;
        a = b;
        b = ctmp;
        #a, b = b, a
        alen, blen = blen, alen
        
    # short circuit.
    if blen - alen >= limit:
        return blen - alen

    cdef char *m1 = <char *>calloc(blen + 1, sizeof(char))
    cdef char *m2 = <char *>malloc((blen + 1) * sizeof(char))
    cdef char *m3 = <char *>malloc((blen + 1) * sizeof(char))
    
    for i from 0 <= i <= blen:
        m2[i] = i 

    for i from 1 <= i <= alen:
        m1[0] = i + 1
        achr = a[i - 1]
        for j from 1 <= j <= blen:
            bchr = b[j- 1]
            #cost = 0 if achr == bchr else 1
            if achr == bchr:
            #if cost == 0:
                m1[j] = m2[j - 1]
            else:
                m1[j] = 1 + imin(m1[j - 1], m2[j - 1], m2[j])

            if i != 1 and j != 1 and achr == b[j - 2] and bchr == a[i - 2]: # and m1[j] > m3[j - 1] + cost:
                    m1[j] = m3[j - 1] #+ cost

        m1, m2 = m2, m1
        strcpy(m3, m2)

    retval = m2[blen]
    free(m3)
    free(m1)
    free(m2)
    return retval

cdef inline int get_arc_dist(arc):
    cdef int arcd
    arcd = <int>PyTuple_GET_ITEM(arc, 1)
    Py_INCREF(arcd)
    return arcd


cdef class BKTree:
    """
    a simple http://en.wikipedia.org/wiki/BK-tree
    see tests for example usage.
    takes a list of words and creates a tree that can be searched
    quickly.  search with tree.find("word", distance)

    follows bearophiles code from here:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/572156

    >>> words = ("abcd", "def", "acdf", "gef", "wersdf", "asdfewd"
    ...        , "sdfwe", "sdfwef", "kljd", "oiouoij", "ououoji", "ouonkoj")

    >>> tree = BKTree(words)
    >>> for i in range(1, 4):
    ...     i, tree.find("ddd", i)
    1 []
    2 ['def']
    3 ['abcd', 'def', 'gef', 'acdf', 'kljd']
    """
    """
    TODO:
    allow a tuple of [(info, word),...]
    >>> words = [('some_name', "word"), 

    cdef object root
    cdef object nodes

    def __init__(self, words):
        root = words[0]
        self.root = root
        self.nodes = {root: []}
        for w in words[1:]:
            if w in self.nodes: continue
            self.addNode(root, w)

    cpdef addNode(self, char *root, char *word):
        cdef int d = edit_distance(root, word, 1000)
        cdef int arc_dist
        if d == 0: return None

        for arc in self.nodes[root]:
            arc_dist = get_arc_dist(arc)

            if d != arc_dist: continue
            self.addNode(arc[0], word)
            break

        else:
            if not word in self.nodes:
                self.nodes[word] = []
            self.nodes[root].append((word, d))


    def find(self, char *word, thresh):
        cdef object results = []
        self._find(self.root, word, thresh, results)
        return results

    cdef void _find(self, char *aroot, char *word, int thresh, results):
        cdef int d = edit_distance(word, aroot, 1000)
        cdef int dmin = d - thresh
        cdef int dmax = d + thresh
        cdef int arc_dist

        if d <= thresh:
            results.append(aroot)

        # todo PyList_GetItem...
        for arc in self.nodes[aroot]:
            arc_dist = get_arc_dist(arc)
            if dmin <= arc_dist <= dmax:
                self._find(arc[0], word, thresh, results)
