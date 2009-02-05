import numpy as np
from npytxt import savetxt, loadtxt, loadgff
from cStringIO import StringIO

import unittest
import os
import time
import sys
import shutil

class NpyTxtTester(unittest.TestCase):
    def setUp(self):
        d = StringIO("M 21 72.1\nF 35 58.33")
        self.a =  np.loadtxt(d, dtype={'names' : ('gender', 'age', 'weight')
                       , 'formats': ('S1', 'i4', 'f8')})
    
        savetxt("g.txt", self.a)

    def test_save(self):
        self.assert_(os.path.exists("g.txt"))
        self.assert_(os.path.exists("g.txt.npy"))

        txt = open("g.txt")
        self.assertEqual(len(txt.readlines()), 3)

    def test_load(self):
        b = loadtxt("g.txt")
        self.assertEqual(self.a.dtype, b.dtype, (self.a.dtype, b.dtype))
        self.assertEqual(self.a[1], b[1])

        # can extract by column
        self.assertEqual(b[b['gender'] == 'M'].shape[0], 1)


    def test_names(self):
        # set the names manually, and it will figure out the format string. 
        # from the dtype of the array.
        savetxt("g.txt", self.a, names=("a", "b", "c"))
        c = loadtxt("g.txt")
        self.assertEqual(["a", "b", "c"], sorted(c.dtype.fields.keys()))

        self.assertEqual(self.a[1]['age'], c[1]['b'])
        self.assertEqual(self.a[1]['gender'], c[1]['a'])
        self.assertEqual(c[0]['c'], 72.1, c[0])


    def test_change(self):
        d = loadtxt("g.txt")
        assert d.shape[0] == 2
        self.assertEqual(d.dtype, self.a.dtype)
        
        # add a column and make sure the npy is upated.
        # os.stat has a resolution of 1 second?
        time.sleep(1)
        shutil.copyfile("g.txt", "h.txt")
        fh = open('h.txt', 'a')
        fh.write("M\t33\t21.99\n")
        fh.close()
        time.sleep(1)
        
        d = loadtxt("h.txt")
        self.assertEqual(d.dtype, self.a.dtype)
        assert d.shape[0] == 3, (d.shape[0], open("h.txt").readlines())
        assert d[0]['gender'] == 'M'

        e = np.load('g.txt.npy')
        assert e[0]['weight'] == d[0]['weight'], (d, e)
        os.unlink("h.txt.npy")
        os.unlink("g.txt.npy")
        os.unlink("h.txt")

    def tearDown(self):
        os.unlink("g.txt")


class NpyGFFTester(unittest.TestCase):
    def setUp(self):
        gffstr = """\
##gff-version 3
1\tucb\tgene\t2234602\t2234702\t.\t-\t.\tID=grape_1_2234602_2234702;match=EVM_prediction_supercontig_1.248,EVM_prediction_supercontig_1.248.mRNA
1\tucb\tgene\t2300292\t2302123\t.\t+\t.\tID=grape_1_2300292_2302123;match=EVM_prediction_supercontig_244.8
1\tucb\tgene\t2303615\t2303967\t.\t+\t.\tID=grape_1_2303615_2303967;match=EVM_prediction_supercontig_244.8
1\tucb\tgene\t2303616\t2303966\t.\t+\t.\tParent=grape_1_2303615_2303967
1\tucb\tgene\t3596400\t3596503\t.\t-\t.\tID=grape_1_3596400_3596503;match=evm.TU.supercontig_167.27
1\tucb\tgene\t3600651\t3600977\t.\t-\t.\tmatch=evm.model.supercontig_1217.1,evm.model.supercontig_1217.1.mRNA
"""
        d = StringIO(gffstr)
        self.a = loadgff(d)
        #print >>sys.stderr, self.a.dtype

    def test_cols(self):
        self.assertEqual(self.a[0]['seqid'], '1')
        self.assert_(np.all(self.a['end'] == np.array([2234702, 2302123, 2303967,
                                                2303966, 3596503, 3600977])))

    def test_attrs(self):
        self.assertEqual(self.a[1]['attrs']['ID'], 'grape_1_2300292_2302123')
        #self.assertEqual(self.a[1]['attrs']['match'], 'EVM_prediction_supercontig_244.8')

def main():
    unittest.main()
    tearDown()

if __name__ == "__main__":
    main()
