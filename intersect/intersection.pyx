#!/usr/bin/python
"""
`Intersecter` data structure for a set of intervals which supports fast 
overlap queries and preserves all information about the intervals (unlike
bitset projection methods). 
"""
import operator

__all__ = ["Feature", "Intersecter", "distance"]

# NOTE: can run docstrings directly from this file using:
# nosetests --with-doctest --doctest-extension=.pyx intersection.pyx

cdef class Feature:
    """\
    Basic feature, with required integer start and stop properties.
    Also accpets optional strand as +1 or -1 (used for up/downstream queries),
    a name, and any arbitrary data is sent in on the info keyword argument

    >>> from intersection import Feature

    >>> f1 = Feature(23, 36)
    >>> f2 = Feature(34, 48, strand=-1, name="fred", info={'chr':12, 'anno':'transposon'})
    >>> f2
    Feature(34, 48, strand=-1, name="fred", {'anno': 'transposon', 'chr': 12})

    """
    cdef public int start, stop, strand
    cdef public object info, name

    def __cinit__(self, int start, int stop, int strand=0, name="", object info=None):
        assert start <= stop, "start must be less than stop"
        self.start  = start
        self.stop   = stop
        self.strand = strand
        self.name   = name
        self.info   = info

    def __repr__(self):
        fstr = "Feature(%d, %d" % (self.start, self.stop)
        if self.strand != 0:
            fstr += ", strand=%d" % self.strand
        if len(self.name):
            fstr += ', name="' + str(self.name) + '"'
        if not self.info is None:
            fstr += ", " + str(self.info)
        fstr += ")"
        return fstr


# like python's bisect_left. find the _lowest_ index where the value x 
#could be  inserted to maintain order in the list intervals
cdef int binsearch_left_start(list intervals, int x, int lo, int hi):
    cdef int mid
    cdef Feature f
    while lo < hi:
        mid = (lo + hi)//2
        f = intervals[mid]
        if f.start < x: lo = mid + 1
        else: hi = mid
    return lo
    
# like python's bisect_right find the _highest_ index where the value x 
# could be  inserted to maintain order in the list intervals
cdef int binsearch_right_stop(list intervals, int x, int lo, int hi):
    cdef int mid
    cdef Feature f
    while lo < hi:
        mid = (lo + hi)/2
        f = intervals[mid]
        if x < f.start: hi = mid
        else: lo = mid + 1
    return lo


# hack since pyrex doesnt support closures.
# used to sort the features by the distance from a query feature.
def _dist_compare(a, b):
    return cmp(distance(a[0], a[1]), distance(b[0], b[1]))


cpdef int distance(Feature f1, Feature f2):
    """\
    Distance between 2 features. The integer result is always positive or zero.
    If the features overlap or touch, it is zero.
    >>> from intersection import Feature, distance
    >>> distance(Feature(1, 2), Feature(12, 13))
    10
    >>> distance(Feature(1, 2), Feature(2, 3))
    0
    >>> distance(Feature(1, 100), Feature(20, 30))
    0

    """
    if f1.stop < f2.start: return f2.start - f1.stop
    if f2.stop < f1.start: return f1.start - f2.stop
    return 0


cdef class Intersecter:
    """\
    Data structure for performing intersect and neighbor queries on a 
    set of intervals. Algorithm uses simple binary search along with
    knowledge of the longest interval to perform efficient queries. 

    Usage
    =====
    >>> from intersection import Intersecter, Feature
    >>> intersecter = Intersecter()

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, name, and info.

    >>> Feature(1, 22, strand=-1, name="fred", info={'chr':12, 'anno': 'anything'})
    Feature(1, 22, strand=-1, name="fred", {'anno': 'anything', 'chr': 12})


    >>> intersecter.add_interval(Feature(0, 10, -1))
    >>> intersecter.add_interval(Feature(3, 7, 1))
    >>> intersecter.add_interval(Feature(3, 40, -1))
    >>> intersecter.add_interval(Feature(13, 50, 1))

    Queries
    -------

    find
    ++++

    >>> intersecter.find(2, 5)
    [Feature(0, 10, strand=-1), Feature(3, 7, strand=1), Feature(3, 40, strand=-1)]
    >>> intersecter.find(11, 100)
    [Feature(3, 40, strand=-1), Feature(13, 50, strand=1)]
    >>> intersecter.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> intersecter.left(Feature(0, 1))
    []
    >>> intersecter.left(Feature(11, 12))
    [Feature(0, 10, strand=-1)]


    up/downstream
    +++++++++++++
    up/downstream method behave exactly like left/right, except that
    the direction is determined by the strand of the query feature. 
    If the strand is 1, then upstream is left, downstream is right.

    If the strand is -1, then upstream is right, downstream is left.
    >>> intersecter.upstream(Feature(11, 12, strand=1))
    [Feature(0, 10, strand=-1)]
    >>> intersecter.upstream(Feature(11, 12, strand=-1))
    [Feature(13, 50, strand=1)]

    all of these method take an argument 'n' for the number of results desired.
    >>> intersecter.upstream(Feature(1, 2, strand=-1), n=3)
    [Feature(3, 7, strand=1), Feature(3, 40, strand=-1), Feature(13, 50, strand=1)]
    
    nearest neighbors
    +++++++++++++++++
    >>> intersecter.nearest_neighbors(Feature(1, 2))
    [Feature(0, 10, strand=-1)]

    >>> intersecter.nearest_neighbors(Feature(1, 2), n=2)
    [Feature(0, 10, strand=-1), Feature(3, 7, strand=1)]

    """

    # since intervals are sorted by start, also have to know the max_len (see find)
    cdef int max_len
    # if an item is added, the list must be resorted. (see # add_interval)
    cdef int dirty
    # a list storing the intervals/Feature's.
    cdef object intervals

    # ---- Basic API --------------------------------------------------

    def __init__(self, intervals=None):
        self.max_len = 0
        self.dirty = 0

        if not intervals is None:
            self.intervals = intervals
            self.intervals.sort(key=operator.attrgetter('start'))
            self.max_len = max([i.stop - i.start for i in intervals])
            if self.max_len < 1: self.max_len = 1
        else:
            self.intervals = []

    def add_interval(self, interval, insort=False): # add Feature type slow things down
        """Add a feature to the intersecter, put it in order by start.
        if insert is True, it will insert the feature into it's sorted
        location in the Intersector. This is only good for inserting
        few values into an existing tree.
        When doing a large number of inserts, always use insort=False.
        """
        cdef int ilen = interval.stop - interval.start
        if ilen > self.max_len: self.max_len = ilen
        assert ilen > -1, ("Features must have length >= 1", interval) 
        if insort:
            idx = binsearch_left_start(self.intervals, interval.start, 0, len(self.intervals))
            self.intervals.insert(idx, interval)
        else:
            self.intervals.append(interval)
            self.dirty = 1

    cdef void prepare(Intersecter self):
        # internal. the search needs the intervals to be sorted by start
        self.intervals.sort(key=operator.attrgetter('start'))
        self.dirty = 0

    def find(Intersecter self, int start, int stop):
        """Return a object of all stored intervals intersecting between (start, end) inclusive."""
        if self.dirty == 1: self.prepare()
        cdef list intervals = self.intervals
        cdef int ilen = len(intervals)
        # NOTE: we only search for starts, since any feature that starts within max_len of
        # the query could overlap, we must subtract max_len from the start to get the needed
        # search space. everything else proceeds like a binary search.
        # (but add distance calc for candidates).
        cdef int ileft  = binsearch_left_start(intervals, start - self.max_len, 0, ilen)
        cdef int iright = binsearch_right_stop(intervals, stop, ileft, ilen)
        cdef Feature query = Feature(start, stop)
        # we have to check the distance to make sure we didnt pick up anything 
        # that started within max_len, but wasnt as long as max_len
        return [f for f in intervals[ileft:iright] if distance(f, query) == 0]

    cpdef list left(self, Feature f, int n=1):
        """return the nearest n features strictly to the left of a Feature f.
        Overlapping features are not considered as to the left.
   
        f: a Feature object
        n: the number of features to return
        """ 
        if self.dirty == 1: self.prepare()
        cdef list intervals = self.intervals
        cdef int iright = binsearch_left_start(intervals, f.start, 0 , len(intervals)) + 1
        cdef int ileft  = binsearch_left_start(intervals, f.start - self.max_len - 1, 0, iright - 1)  
        cdef Feature other

        results = [(other, f) for other in intervals[ileft:iright] if other.stop < f.start and distance(f, other) != 0]
        results.sort(cmp=_dist_compare)
        results = [r[0] for r in results]
        if len(results) == n: return results

        # have to do some extra work here since intervals are sorted
        # by starts, and we dont know which stops may be around...
        # in this case, we got some extras, just return as many as
        # needed once we see a gap in distances.
        for i in range(n, len(results)):
            if distance(f, results[i - 1]) != distance(f,  results[i]):
                return results[:i]

        # here made it all the way to zero.
        if ileft == 0: return results

        # here, didn't get enough, so move left and try again. 
        # TODO: add tests for this case..
        #return results + self.left(Feature(start=f.start - self.max_len - 2, stop=f.stop), n=n - len(results))
        if len(results) != 0:
            return results + self.left(Feature(start=results[-1].start, stop=f.stop), n=n - len(results))
        return results + self.left(Feature(start=f.start - 1, stop=f.stop), n=n - len(results))


    cpdef list right(self, Feature f, int n=1):
        """return the nearest n features strictly to the right of a Feature f.
        Overlapping features are not considered as to the right.
   
        f: a Feature object
        n: the number of features to return
        """ 
        if self.dirty == 1: self.prepare()
        cdef list intervals = self.intervals
        cdef int ilen = len(intervals)
        cdef int iright = binsearch_right_stop(intervals, f.stop, 0, ilen)
        cdef Feature other
        cdef list results = []
        cdef int i


        while iright < ilen:
            i = len(results)
            if i > n: 
                if distance(f, results[i - 1]) != distance(f, results[i - 2]): 
                    return results[:i - 1]
            other = intervals[iright]
            iright += 1
            if distance(other, f) == 0: continue
            results.append(other)
        return results


    def upstream(self, Feature f, int n=1):
        """find n upstream features where upstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        """
        if f.strand == -1:
            return self.right(f, n)
        return self.left(f, n)


    def downstream(self, Feature f, int n=1):
        """find n downstream features where downstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        """
        if f.strand == -1:
            return self.left(f, n)
        return self.right(f, n)

    def nearest_neighbors(self, Feature f, int n=1, inclusive=True):
        """return the n nearest neighbors to the given feature
        f: a Feature object
        n: the number of features to return
        inclusive: True/False: other operations search only _around_ the query feature,
            if inclusive is False, this query will not find a feature that is completely
            contained within the query feature (this would be desired for up/down stream, 
            but less intuitive for a neighborhood search

            can nearly always ignore this and leave as True
        """
        if self.dirty == 1: self.prepare()
        cdef list intervals = self.intervals
        cdef int ilen = len(intervals)
        cdef int ileft  = binsearch_left_start(intervals, f.start - self.max_len, 0, ilen)
        cdef int iright = binsearch_right_stop(intervals, f.stop, ileft, ilen)
        # TODO: if both left and right dists are gt than the last,
        # then return...
        # TODO: add tessts for overlaps.
        while ileft > 0 or iright < ilen:
            if iright - ileft >= n and distance(f, intervals[iright]) != distance(f, intervals[ileft]):
                return [ff[1] for ff in sorted([(f, iv) for iv in intervals[ileft:iright]], cmp=_dist_compare)]
            ileft -= 1
            iright += 1
            if ileft == -1: ileft = 0
            if iright > ilen: iright = ilen
        return [ff[1] for ff in sorted([(f, iv) for iv in intervals[ileft:iright]], cmp=_dist_compare)]

