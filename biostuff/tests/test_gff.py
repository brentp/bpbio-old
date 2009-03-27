from biostuff.gff_reader import GFFLine, GFFNode
import sys


def test_gff_node():
    f = "tests/data/th.gff" 

    node = GFFNode.yield_nodes(f)
    print >>sys.stderr, node.next(), "\n"
    print >>sys.stderr, node.next(), "\n"
    print >>sys.stderr, node.next(), "\n"
    print >>sys.stderr, node.next(), "\n"
    print >>sys.stderr, node.next(), "\n"


    #yield check_type, blasts, ('eval', 'score', 'pctid'), float
    #yield check_type, blasts, ('query', 'subject'), str
    

def check_type(blasts, attrs, klass):
    for b in blasts:
        for attr in attrs:
            assert isinstance(getattr(b, attr), klass)

