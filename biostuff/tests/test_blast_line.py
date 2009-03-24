from biostuff import BlastLine


def test_blastline():
    f = "tests/data/tabd.blast" 
    blasts = []
    for line in open(f):
        blasts.append(BlastLine(line))

    yield check_type, blasts, ('qstart', 'qstop', 'sstart', 'sstop',
        'nmismatch', 'ngaps'), int

    yield check_type, blasts, ('eval', 'score', 'pctid'), float
    yield check_type, blasts, ('query', 'subject'), str
    

def check_type(blasts, attrs, klass):
    for b in blasts:
        for attr in attrs:
            assert isinstance(getattr(b, attr), klass)

