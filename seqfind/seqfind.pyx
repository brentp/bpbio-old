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



    cdef int alen = strlen(a), blen = strlen(b), retval
    cdef char *ctmp
    cdef size_t i, j
    cdef size_t achr, bchr
     
    if strcmp(a, b) == 0:
        return 0

    if alen > blen:
        ctmp = a;
        a = b;
        b = ctmp;
        alen, blen = blen, alen
        
    # short circuit.
    if blen - alen >= limit:
        return blen - alen

    cdef char *m1 = <char *>calloc(blen + 2, sizeof(char))
    cdef char *m2 = <char *>calloc(blen + 2, sizeof(char))
    cdef char *m3 = <char *>malloc((blen + 2) * sizeof(char))
    
    for i from 0 <= i <= blen:
        m2[i] = i 

    for i from 1 <= i <= alen:
        m1[0] = i + 1
        achr = a[i - 1]
        for j from 1 <= j <= blen:
            bchr = b[j- 1]
            if achr == bchr:
                m1[j] = m2[j - 1]
            else:
                m1[j] = 1 + imin(m1[j - 1], m2[j - 1], m2[j])

            if i != 1 and j != 1 and achr == b[j - 2] and bchr == a[i - 2]: # and m1[j] > m3[j - 1] + cost:
                    m1[j] = m3[j - 1] 

        m1, m2 = m2, m1
        strcpy(m3, m2)

    retval = <int>m2[blen]
    free(m3)
    free(m1)
    free(m2)
    return retval


cdef class Word:
    cdef public char* word
    cdef public object info
    cdef public int distance

    def __cinit__(self, char * word, object info=None): #, char* info):
        self.word = word
        self.info = info

    def __repr__(self):
        return "<Word. {'word': '%s' }>" % self.word

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

    the words can also be a Word() object. which takes the word to
    compare, and any python object:

        >>> words = [Word("asdf", {'info':"this is some info})
        ...            , Word("asfd", {'zombies': [24, 'g']})]
        >>> tree = BKTree(words)

    this is useful for storing any info that goes with the object.
        >>> f = tree.find("asdf, 0)
        >>> f.info['info']
        'zombies'
        >>> f.info['zombies']
        [24, 'g']
    """

    cdef Word root
    cdef public object nodes

    def __init__(self, words):
        cdef int is_word = 0
        cdef Word root 

        if isinstance(words[0], Word):
            is_word = 1
            root = words[0]
        else:
            root = Word(words[0])
        

        cdef Word aword
        self.root = root
        self.nodes = {root: []}
        for w in words[1:]:
            if is_word == 0:
                aword = Word(w)
            else: 
                aword = w
            if aword in self.nodes: continue
            self.addNode(root, aword)

    cdef void addNode(BKTree self, Word root, Word word):
        cdef int d = edit_distance(root.word, word.word, 1000)
        cdef Word arc

        #if d != 0: 
        for arc in self.nodes[root]:
            if d != arc.distance: continue
            self.addNode(arc, word)
            break

        else:
            if not word in self.nodes:
                self.nodes[word] = []
            word.distance = d
            self.nodes[root].append(word)


    cpdef find(self, char *word, int thresh):
        cdef list results = []
        self._find(self.root, Word(word), thresh, results)
        return results

    cdef void _find(BKTree self, Word aroot, Word word, int thresh, list results):
        cdef int arc_dist, d = edit_distance(word.word, aroot.word, 1000)
        cdef int dmin = d - thresh
        cdef int dmax = d + thresh
        cdef Word arc

        if d <= thresh:
            results.append(aroot)

        # todo PyList_GetItem...
        for arc in self.nodes[aroot]:
            arc_dist = arc.distance
            if dmin <= arc_dist <= dmax:
                self._find(arc, word, thresh, results)

