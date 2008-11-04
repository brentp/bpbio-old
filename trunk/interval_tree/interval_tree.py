class IntervalTree(object):
    __slots__ = ('intervals', 'left', 'right', 'center')

    def __init__(self, intervals, depth=16, minbucket=128, _extent=None):
        """\
        `intervals` a list of intervals *with start and stop* attributes.
        `depth` the depth of the tree
        `minbucket` if any node in the tree has fewer than minbucket
                    elements, make it a leaf node

        depth and minbucket usually do not need to be changed. if
        dealing with large numbers (> 1M) of intervals, the depth could
        be increased to 24.

        Useage:

         >>> ivals = [Interval(2, 3), Interval(1, 8), Interval(3, 6)]
         >>> tree = IntervalTree(ivals)
         >>> tree.find(1, 2)
         [Interval(2, 3), Interval(1, 8)]

        this provides an extreme and satisfying performance improvement
        over searching manually over all 3 elements in the list (like
        a sucker). 

        NOTE: any object with start and stop attributes can be used
        in the incoming intervals list.
        """ 

        depth -= 1
        if depth == 0 or len(intervals) < minbucket:
            self.intervals = intervals
            self.left = self.right = None
            return 

        left, right = _extent or \
               (min(i.start for i in intervals), max(i.stop for i in intervals))
        center = (left + right) / 2.0

        
        self.intervals = []
        lefts, rights  = [], []
        

        for interval in intervals:
            if interval.stop <= center:
                lefts.append(interval)
            elif interval.start >= center:
                rights.append(interval)
            else: # overlapping.
                self.intervals.append(interval)
                
        self.left   = lefts  and IntervalTree(lefts,  depth, minbucket, (left,  center)) or None
        self.right  = rights and IntervalTree(rights, depth, minbucket, (center, right)) or None
        self.center = center
 
 
    def find(self, start, stop):
        """find all elements between (or overlapping) start and stop"""
        overlapping = [i for i in self.intervals if i.stop >= start 
                                              and i.start <= stop]

        if self.left and start <= self.center:
            overlapping += self.left.find(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find(start, stop)

        return overlapping

class Interval(object):
    def __init__(self, start, stop):
        self.start = start
        self.stop  = stop
    def __repr__(self):
        return "Interval(%i, %i)" % (self.start, self.stop)

if __name__ == '__main__':
    

    def brute_force_find(intervals, start, stop):
        return [i for i in intervals if i.stop >= start and i.start <= stop]


    import random, time
    def rand():
        s = random.randint(1, 4000000)
        return Interval(s, s + random.randint(2000, 6000))
    intervals = [rand() for i in xrange(30000)]
    START, STOP = 300000, 400000
    intervals.append(Interval(0, 500000))

    tree = IntervalTree(intervals)
    t = time.time()
    res = tree.find(START, STOP)
    treetime = time.time() - t
    t = time.time()
    bf = [i for i in intervals if i.stop >= START and i.start <= STOP]
    btime = time.time() - t
    assert len(bf) == len(res) , (len(bf), len(res), set(bf).difference(res), START, STOP)
    print treetime, btime, btime/treetime

    import doctest
    doctest.testmod()

