from biostuff import BlastLine

some_attrs = ('qstart', 'qstop', 'sstart', 'sstop', 'pctid', 'score', 'query',
        'subject')

def test_blastline():
    f = "tests/data/tabd.blast" 
    blasts = []
    for line in open(f):
        blasts.append(BlastLine(line))

    yield check_type, blasts, ('qstart', 'qstop', 'sstart', 'sstop',
        'nmismatch', 'ngaps'), int

    yield check_type, blasts, ('evalue', 'score', 'pctid'), float
    yield check_type, blasts, ('query', 'subject'), str
    

def check_type(blasts, attrs, klass):
    for b in blasts:
        for attr in attrs:
            assert isinstance(getattr(b, attr), klass)


def test_to_string():
    f = "tests/data/tabd.blast" 
    for line in open(f):
        a = BlastLine(line)
        b = BlastLine(a.to_blast_line())

        # works better than string comparison because of floats.
        for attr in some_attrs:
            assert getattr(a, attr) == getattr(b, attr), (a, b, attr)
