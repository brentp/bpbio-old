from pyfasta import Fasta
import unittest

class FastaTest(unittest.TestCase):

    def setUp(self):
        self.f = Fasta('tests/data/three_chrs.fasta')

    def test_keys(self):

        self.assertEqual(sorted(self.f.keys()), ['chr1', 'chr2', 'chr3'])

    
    def test_mmap(self):
        seq = self.f['chr2']
        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[1], 'A')

        self.assertEqual(seq[-1], 'T')
        self.assertEqual(seq[-2], 'A')

        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[6:9], 'AAA')

    def test_shape(self):
        self.assertEqual(len(self.f['chr2']), 80)
        self.assertEqual(len(self.f['chr3']), 3600)

    def test_tostring(self):
        s = 'TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT'
        self.assertEqual(str(self.f['chr2']), s)


    def test_slice(self):
        f = self.f
        self.assertEqual(str(f['chr3'][:5]), 'ACGCA')

    
        seq = 'TACGCACGCTAC'
        self.assertEqual(seq, f['chr3'][-12:])

        self.assertEqual(f['chr3'][0:5][::-1], 'ACGCA')

    def tearDown(self):
        import os
        os.unlink('tests/data/three_chrs.fasta.gz')
        os.unlink('tests/data/three_chrs.fasta.gdx')


if __name__ == "__main__":
    unittest.main()
