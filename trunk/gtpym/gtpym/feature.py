"""
# optional gff file as first arg.
>>> fi = FeatureIndexMemory('tests/data/grape.gff3')
>>> for i, (fid, feature) in enumerate(fi.iteritems()):
...    if i > 3: break
...    feature
FeatureNode(start=163062, end=163582, seqid="1")
FeatureNode(start=287391, end=289450, seqid="1")
FeatureNode(start=26433, end=45352, seqid="1")
FeatureNode(start=101683, end=103732, seqid="1")

>>> fi.seqids
['1']


# by feature.attribs["ID"]
>>> f = fi['GSVIV00029993001']
>>> f
FeatureNode(start=4764, end=20015, seqid="1")

>>> FeatureIndexMemory.by_types(f)
{'UTR': [FeatureNode(start=4764, end=4978, seqid="1"), FeatureNode(start=6832, end=6870, seqid="1")], 'gene': [FeatureNode(start=4764, end=20015, seqid="1")], 'mRNA': [FeatureNode(start=4764, end=20015, seqid="1")], 'CDS': [FeatureNode(start=6871, end=6893, seqid="1"), FeatureNode(start=7031, end=7153, seqid="1"), FeatureNode(start=19964, end=20015, seqid="1")]}

#wrapper for get_features_for_range
>>> a = fi[369:3700:'1'][0]
>>> b = fi.get_features_for_range(369, 3700, '1')[0]
>>> a.start == b.start and a.end == b.end and a.get_attribute("ID") == b.get_attribute("ID")
True


>>> fi[369:3700:'1']
[FeatureNode(start=1579, end=2025, seqid="1")]


# by seqid -- wrapper for get_features_for_seqid
>>> a = fi['1'][:2]
>>> b = fi.get_features_for_seqid("1")[:2]
>>> [f.start for f in a] == [f.start for f in b]
True

>>> FeatureNode('6', 'gene', 123, 234, '+')
FeatureNode(start=123, end=234, seqid="6")

"""
import gt
import collections
import os

__all__ = ['FeatureIndexMemory', 'FeatureNode']

class FeatureNode(gt.FeatureNode):

    def __new__(cls, seqid, type, start, end, strand):
        obj = gt.FeatureNode.create_new(seqid, type, start, end, strand)
        #obj.attribs = dict(obj.each_attribute())
        return obj

    @classmethod
    def create_from_ptr(self, ptr, newref=False):
        fn = gt.FeatureNode.create_from_ptr(ptr)
        #fn.attribs = dict(fn.each_attribute())
        return fn


class FeatureIndexMemory(gt.FeatureIndexMemory):
    def __init__(self, *gff3files, **kwargs):
        with_ids = kwargs.get('with_ids', True)
        gt.FeatureIndexMemory.__init__(self)
        self.seqids = set([])
        for gff3file in gff3files:
            assert os.path.exists(gff3file), (gff3file, "does not exist")
            self.add_gff3file(gff3file)
        self.seqids = sorted(self.get_seqids())
        self._ids = {}
        if with_ids:
            self._load_ids()

    @staticmethod
    def by_types(f):
        types = collections.defaultdict(list)
        for feature_bit in f(method="depth_first"):
            types[feature_bit.type].append(feature_bit)
        return dict(types)

    def __getitem__(self, fid):
        if isinstance(fid, slice):
            return self.get_features_for_range(fid.start, fid.stop, str(fid.step))
        else:
            try:
                ptr = self._ids[fid]
                return FeatureNode.create_from_ptr(ptr, newref=True)
            except KeyError:
                if fid in self.seqids:
                    return self.get_features_for_seqid(fid)
                else:
                    raise

    def _load_ids(self):
        self._ids = {}
        for seqid in self.seqids:
            for f in self.get_features_for_seqid(seqid):
                self._ids[f.get_attribute('ID')] = f._as_parameter_

    def iteritems(self):
        for fid, ptr in self._ids.iteritems():
            yield fid, FeatureNode.create_from_ptr(ptr)
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
