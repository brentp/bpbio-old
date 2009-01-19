import numpy as np
from ncsv import savetxt, loadtxt
from cStringIO import StringIO

import unittest
import os

class NCSVTester(unittest.TestCase):
    def setUp(self):
        d = StringIO("M 21 72.1\nF 35 58.33")
        self.a =  np.loadtxt(d, dtype={'names' : ('gender', 'age', 'weight')
                       , 'formats': ('S1', 'i4', 'f8')})
    
        savetxt(self.a, "g.txt")

    def test_save(self):
        self.assert_(os.path.exists("g.txt"))

        txt = open("g.txt")
        self.assertEqual(len(txt.readlines()), 3)

    def test_load(self):
    
        b = loadtxt("g.txt")
        self.assertEqual(self.a.dtype, b.dtype)
        self.assertEqual(self.a[1], b[1])


    def test_names(self):
        # set the names manually, and it will figure out the format string. 
        # from the dtype of the array.
        savetxt(self.a, "g.txt", names=("a", "b", "c"))
        c = loadtxt("g.txt")
        self.assertEqual(["a", "b", "c"], sorted(c.dtype.fields.keys()))

        self.assertEqual(self.a[1]['age'], c[1]['b'])
        self.assertEqual(c[0]['c'], 72.1, c[0])

    def tearDown(self):
        os.unlink("g.txt")

def main():
    unittest.main()

if __name__ == "__main__":
    main()
