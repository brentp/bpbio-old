from psyco.classes import *
from misc import read_gff

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from intersection import Feature
from intersection import Intersecter as CIntersecter
from bx.intervals.intersection import Intersecter
import time
import random

TEST = len(sys.argv) > 1
class PyFeature:
    def __init__(self, start, stop, *args):
        self.start = start
        self.stop  = stop
        self.end   = stop
        self.stuff = args
    def __repr__( self ):
        return "%s( %d, %d, %s )" % (self.__class__.__name__,  self.start, self.stop, repr( self.stuff ) )


fname = "tests/rice_condensed.gff"
all_chrs = read_gff(fname, name_fn=lambda line: line[8].split(" ")[1])
data = all_chrs['1']


rands = [random.randint(0, 10000000) for i in range(50000)]

inter = Intersecter()
cinter = CIntersecter()


t0 = time.time()
for feat in data:
    inter.add_interval(PyFeature(feat['start'], feat['stop'], feat['strand'], feat['name'], feat))

t1 = time.time()

for feat in data:
    cinter.add_interval(Feature(feat['start'], feat['stop'], feat['strand'], feat['name'], feat))

t2 = time.time()
print "time to create (python)", t1 - t0
print "time to create (cython)", t2 - t1
print "speedup", (t1 - t0)/(t2 - t1)
print "...searching an intersecter with %i features %i times..." % (len(data), len(rands))
#cinter.find(1, 2)
#inter.find(1,2)

# add a feature that's the same extent as all the features.
# this simulates the worst-case scenario for the Cython version.
#cinter.add_interval(Feature(1, 10000000/1.))
#inter.add_interval(PyFeature(1, 10000000/1.))
t2 = time.time()


size = 5000

pfound = []
for r in rands:
    found = inter.find(r, r + size)
    if TEST: pfound.append(([(f.start, f.stop) for f in found], (r, r + size)))

t3 = time.time()


cfound = []
for r in rands:
    found = cinter.find(r, r + size)
    if TEST: cfound.append(([(f.start, f.stop) for f in found], (r, r + size)))

t4 = time.time()

print
print "time to search (python)", t3 - t2
print "time to search (cython)", t4 - t3
print "speedup", (t3 - t2)/(t4 - t3)

if TEST:
    print 'pfound[0]',sorted(pfound[0])
    print 'cfound[0]',sorted(cfound[0])

for p, c in zip(pfound, cfound):
    c0 = sorted(c[0])
    p0 = sorted(p[0])
    if len(c0) != len(p0): 
        # should only happedn when either the start or end 
        # of the query == start of subject...
        # this is because bx does not consider Feature(23, 45) 
        # and Feature(45, 55) to overlap. i do...

        query = c[1]
        assert c[1] == p[1] 
        diff = list(set(p0).symmetric_difference(c0))[0]
        if diff[0] in query or diff[1] in query: continue
        print 'query:',c[1]
        print 'cython:', c0
        print 'python:', p0
        raise Exception()
