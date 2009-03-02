import os
import cPickle
import mmap # 4G limit?
import string
import gc

_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
complement  = lambda s: s.translate(_complement)

class Fasta(dict):
    def __init__(self, fasta_name, flatten_inplace=True):
        self.fasta_name, self.index = Fasta.flatten(fasta_name, flatten_inplace=flatten_inplace)
        self.fasta_handle = open(self.fasta_name, 'rb+')
        self._chrs = {}
        self.chr = None
        self._load_mmap()

    @classmethod
    def is_up_to_date(klass, a, b):
        # need the + 10 because the index is created after the index file
        return os.path.exists(a) and os.stat(a).st_mtime + 10 > os.stat(b).st_mtime


    @classmethod
    def flatten(self, fasta, flatten_inplace=True):
        """remove all newlines from the sequence in a fasta file"""
        # if save_to is None: modify the file inplace 

        gdx = None
        if not flatten_inplace:
            idx = fasta.rfind('.fa')
            if idx < 0: idx = len(fasta)
            save_to = fasta[:idx] + '.flat' + fasta[idx:]
            gdx = save_to + ".gdx"
            # if the index isnt up to date, reflatten the file.
            if Fasta.is_up_to_date(gdx, fasta):
                return save_to, Fasta._load_index(save_to)


        if flatten_inplace:
            gdx = fasta + '.gdx' 
            # it was modified inplace already.
            if Fasta.is_up_to_date(gdx, fasta):
                return fasta, Fasta._load_index(fasta)
            
            # in order to flatten 'inplace', actually
            # make a tempfile and then overwrite the original.
            import shutil, tempfile
            save_to = tempfile.mktemp('.fasta')

        size = os.path.getsize(fasta)
        fh = open(fasta, 'r+')
        out = open(save_to, 'w')

        mm = mmap.mmap(fh.fileno(), size)

        sheader = mm.find('>')
        snewline = mm.find('\n', sheader)
        idx = {}
        while sheader < len(mm):
            header = mm[sheader:snewline + 1]
            out.write(header)

            sheader = mm.find('>', snewline)
            if sheader == -1: sheader = len(mm)

            seq  = mm[snewline + 1: sheader].replace('\n','')
            out.write(seq + '\n')
            out.flush()
            pos = out.tell()
            idx[header[1:].strip()] = (pos - len(seq) - 1, pos)

            snewline = mm.find('\n', sheader)

        p = open(gdx, 'wb') 
        cPickle.dump(idx, p)
        p.close(); fh.close(); out.close()

        try:
            shutil.copyfile(save_to, fasta)
            os.unlink(save_to)
            save_to = fasta
        except: 
            # it wasn't a tempfile. ignore
            pass

        return save_to, idx

    def _load_mmap(self):
        for name, (start, stop) in self.index.iteritems():
            # lazy, so the memmap isnt loaded until it's requested.
            # TODO. simplify this by just using self.index
            self._chrs[name] = {'offset': start, 'shape': (stop - start -1,) }

    def iterkeys():
        for k in self.keys(): yield k
    def keys(self):
        return self.index.keys()
    
    def __contains__(self, key):
        return key in self.index.keys()

    def __getitem__(self, i):
        # this implements the lazy loading an only allows a single 
        # memmmap to be open at one time.
        assert i in self.index.keys()
        if i != self.chr:
            if self.chr: self.chr.close()
            c = self._chrs[i]
            # http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/core/memmap.py#L18
            #self.chr[i] = numpy.memmap(self.fasta_handle, offset=c['offset'], shape=c['shape'], dtype='c')

            self.chr = mmap.mmap(self.fasta_handle, c['shape'], None,
                    mmap.ACCESS_COPY, c['offset'])
        return self.chr

    def close(self):
        if self.chr: self.chr.close()
        self.fasta_handle.close()

    def __del__(self):
        self.close()

    @classmethod
    def _load_index(self, fasta):
        """ """
        gdx = open(fasta + '.gdx', 'rb')
        return cPickle.load(gdx)

    def sequence(self, f, auto_rc=True
            , exon_keys=None):
        """
        take a feature and use the start/stop or exon_keys to return
        the sequence from the assocatied fasta file:
        f: a feature
        auto_rc : if True and the strand of the feature == -1, return
                  the reverse complement of the sequence

            >>> from genedex import Fasta
            >>> f = Fasta('tests/data/three_chrs.fasta')
            >>> f.sequence({'start':1, 'stop':2, 'strand':1, 'chr': 'chr1'})
            'AC'

        NOTE: these 2 are reverse-complement-ary because of strand
            >>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr1'})
            'CAG'
            >>> f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr1'})
            'CTG'
        
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
            sequence = fasta[f['start'] -1: f['stop']].tostring()

        if auto_rc and f.get('strand') in ('-1', -1, '-'):
            sequence = complement(sequence)[::-1]
        return sequence

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
                seq += fasta[start -1:stop].tostring()
            return seq
        return None

    def iteritems(self):
        for k in self.keys():
            yield k, self[k]



if __name__ == "__main__":
    import doctest
    doctest.testmod()

