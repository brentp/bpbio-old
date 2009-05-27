import os
import cPickle
import numpy as np
import string
import gc
import gzip
import mmap
import sys

_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
complement  = lambda s: s.translate(_complement)


class FastaGz(gzip.GzipFile):
    __slots__ = ('gz_name', 'start', 'stop', 'mode', 'seek', 'read')

    def __init__(self, gz_name, start, stop):
        self.gz_name = gz_name
        self.stop    = stop
        self.start   = start
        gzip.GzipFile.__init__(self, gz_name, mode='rb')

    def __len__(self):
        return self.stop - self.start

    def __getitem__(self, islice):

        if isinstance(islice, int):
            if islice < 0:
                self.seek(self.stop + islice)
            else:
                self.seek(self.start + islice)
            return self.read(1)

        if islice.start == 0 and islice.stop == sys.maxint:
            self.seek(self.start)
            return self.read(self.stop - self.start)

        #if islice.start > islice.stop:
            #    self.seek(self._offset + islice.stop)
            #l = islice.start - islice.stop
            #return complement(self.read(l))[::-1]

        self.seek(self.start + islice.start)
        l = islice.stop - islice.start
        if islice.step in (1, None):
            return self.read(l)

        #elif islice.step == -1:
            #seq = self.read(l)
            #return complement(seq)[::-1]

    def __str__(self):
        return self[:]

    def __repr__(self):
        return "%s('%s', %i..%i)" % (self.__class__.__name__, self.gz_name,
                                   self.start, self.stop)

class Fasta(dict):
    def __init__(self, fasta_name):
        self.fasta_name = fasta_name
        self.gdx = fasta_name + ".gdx"
        self.gz = fasta_name + ".gz"
        self.index = self.gzify()

        self.chr = {}

    @classmethod
    def is_up_to_date(klass, a, b):
        return os.path.exists(a) and os.stat(a).st_mtime > os.stat(b).st_mtime


    def gzify(self):
        """remove all newlines from the sequence in a fasta file"""

        if Fasta.is_up_to_date(self.gdx, self.fasta_name) \
                     and Fasta.is_up_to_date(self.gz, self.fasta_name):
            return Fasta._load_index(self.gdx)


        out = gzip.open(self.gz, 'wb')

        fh = open(self.fasta_name, 'r+')
        mm = mmap.mmap(fh.fileno(), os.path.getsize(self.fasta_name))

        # do the flattening (remove newlines)
        sheader = mm.find('>')
        snewline = mm.find('\n', sheader)
        idx = {}
        pos = 0
        while sheader < len(mm):
            header = mm[sheader:snewline + 1]
            #out.write(header)

            sheader = mm.find('>', snewline)
            if sheader == -1: sheader = len(mm)

            seq  = mm[snewline + 1: sheader].replace('\n','')
            out.write(seq)
            out.flush()
            p1 = out.tell()
            idx[header[1:].strip()] = (pos, p1)
            pos = p1 

            snewline = mm.find('\n', sheader)

        p = open(self.gdx, 'wb')
        cPickle.dump(idx, p)
        p.close(); fh.close(); out.close()
        return idx


    def iterkeys():
        for k in self.keys(): yield k
    def keys(self):
        return self.index.keys()
    
    def __contains__(self, key):
        return key in self.index.keys()

    def __getitem__(self, i):
        # this implements the lazy loading an only allows a single 
        # memmmap to be open at one time.
        if i in self.chr:
            return self.chr[i]

        c = self.index[i]
        self.chr[i] = FastaGz(self.gz, c[0], c[1])
        return self.chr[i]

    @classmethod
    def _load_index(self, path):
        """ """
        gdx = open(path, 'rb')
        try:
            return cPickle.load(gdx)
        finally:
            gdx.close()

    def sequence(self, f, asstring=True, auto_rc=True
            , exon_keys=None):
        """
        take a feature and use the start/stop or exon_keys to return
        the sequence from the assocatied fasta file:
        f: a feature
        asstring: if true, return the sequence as a string
                : if false, return as a numpy array
        auto_rc : if True and the strand of the feature == -1, return
                  the reverse complement of the sequence

            >>> from pyfasta import Fasta
            >>> f = Fasta('tests/data/three_chrs.fasta')
            >>> f.sequence({'start':1, 'stop':2, 'strand':1, 'chr': 'chr1'})
            'AC'

            >>> f.sequence({'start':1, 'stop':2, 'strand': -1, 'chr': 'chr1'})
            'GT'

            >>> f.index
            {'chr3': (160, 3760), 'chr2': (80, 160), 'chr1': (0, 80)}

        NOTE: these 2 are reverse-complement-ary because of strand
        #>>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr1'})
            'CAG'
            >>> f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr1'})
            'CTG'


            >>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr3'})
            'TGC'
            >>> f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr3'})
            'GCA'

            >>> f['chr3'][:][-10:]
            'CGCACGCTAC'

        
        a feature can have exons:
            >>> feat = dict(start=9, stop=19, strand=1, chr='chr1'
            ...    , exons=[(9,11), (13, 15), (17, 19)])

        by default, it just returns the full sequence between start
        and stop.
            >>> f.sequence(feat)
            'ACTGACTGACT'

        but if exon_keys is set to an iterable, it will search for
        those keys and will use the first to create a sequence and
        return the concatenated result.
            >>> f.sequence(feat, exon_keys=('rnas', 'exons'))
            'ACTACTACT'

        Note that sequence is 2 characters shorter than the entire
        feature, to account for the introns at base-pairs 12 and 16.

        Also note, it looks for an item with key of 'rnas', and didn't
        fine one, so it continued on to 'exons'. If it doesn't find
        any of the exon keys, it will fall back on the start, stop of
        the feature:
            >>> f.sequence(feat, exon_keys=('fake', 'also_fake'))
            'ACTGACTGACT'
        """
        assert 'chr' in f and f['chr'] in self, (f, f['chr'], self.keys())
        fasta    = self[f['chr']]
        sequence = None
        if not exon_keys is None:
            sequence = self._seq_from_keys(f, fasta, exon_keys)

        if sequence is None:
            sequence = fasta[(f['start'] - 1): f['stop']]

        if auto_rc and f.get('strand') in (-1, '-1', '-'):
            sequence = complement(sequence)[::-1]

        if asstring: return sequence
        return numpy.array(sequence, dtype='c')

    def _seq_from_keys(self, f, fasta, exon_keys, base='locations'):
        """Internal:
        f: a feature dict
        fasta: a Fasta object
        exon_keys: an iterable of keys, to look for start/stop
                   arrays to get sequence.
        base: if base ('locations') exists, look there fore the
        exon_keys, not in the base of the object:
            {'name': 'OS11G42690', 'stop': 25210251, 'locations':
            {'CDS': [(25210018, 25210251)]}, 'start': 25210018, 'chr':
            '11', 'strand': -1} set(['TRNA', 'CDS'])
        """
        fbase = f.get(base, f)
        for ek in exon_keys:
            if not ek in fbase: continue
            locs = fbase[ek]
            seq = ""
            for start, stop in locs:
                seq += fasta[start -1:stop]
            return seq
        return None

    def iteritems(self):
        for k in self.keys():
            yield k, self[k]



if __name__ == "__main__":
    import doctest
    doctest.testmod()

