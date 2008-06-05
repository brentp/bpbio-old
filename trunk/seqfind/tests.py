#
from seqfind import edit_distance
from seqfind import BKTree, Word
import time
import sys


def test_edit_distance():
    print "testing..."
    assert edit_distance("ab", "ba", 20) == 1

    tests = [["", ""], ["a", ""], ["", "a"], ["a", "a"], ["x", "a"],
              ["aa", ""], ["", "aa"], ["aa", "aa"], ["ax", "aa"], ["a", "aa"], ["aa", "a"],
              ["abcdef", ""], ["", "abcdef"], ["abcdef", "abcdef"],
              ["vintner", "writers"], ["vintners", "writers"]];
    expected = [0, 1, 1, 0, 1, 2, 2, 0, 1, 1, 1, 6, 6, 0, 5, 4]
    for i, (a, b) in enumerate(tests):
        assert edit_distance(a, b, 100) == expected[i], (a, b, expected[i], edit_distance(a, b, 100))
        assert edit_distance(a, b, 100) == edit_distance(b, a, 100)

    transpositions = [("OT", "OST", 1), ("TO", "OST", 3), ("ab", "ba", 1), ("abc", "bac", 1)
            , ("surgery", "survey", 2), ("asdf", "adsf", 1), ("asxxdf", "adxxsf", 2), ("aaaaaaaab", "aaaaaaaba", 1)]
    for a, b, expected in transpositions:
        assert edit_distance(a, b, 100) == expected, (a, b, edit_distance(a, b, 100), expected)
        assert edit_distance(a, b, 100) == edit_distance(b, a, 100)
    print 'done'

    import random
    import string
    word = "".join([random.choice(string.letters) for x in range(10)])
    w2 = word[:][:-2]
    assert edit_distance(word, w2, 10) == 2
    w2 = list(word[:])

    w2[4], w2[5] = w2[5], w2[4]
    w2 = "".join(w2)
    assert edit_distance(word, w2, 100) < 2

    print edit_distance("bzzcbczacd", "zzzzzzcacd", 100)

def test_same_words():
    words = ("abcd", "def", "acdf", "gef", "wersdf", "asdfewd"
            ,"abcd", "def", "acdf", "gef", "wersdf", "asdfewd")
    t = BKTree(words)
    abcds = [w.word for w in t.find("abcd", 0)] 
    assert abcds[0] == abcds[1] == "abcd", abcds
    more = [w.word for w in t.find("abcd", 2)] 
    assert more.count("abcd") == 2, more
    assert more.count("acdf") == 2, more
    


def test_tree():
    words = ("abcd", "def", "acdf", "gef", "wersdf", "asdfewd"
            , "sdfwe", "sdfwef", "kljd", "oiouoij", "ououoji", "ouonkoj")

    tree = BKTree(words)
    print "made tree"
    print tree.nodes

    for i in range(6):
        print i, "ddd", tree.find("ddd", i)

def test_words():
    words = ("abcd", "def", "acdf", "gef", "wersdf", "asdfewd"
            , "sdfwe", "sdfwef", "kljd", "oiouoij", "ououoji", "ouonkoj")

    words = [Word(w, {'something':"asdf", "asdf": [22, 33]}) for w in words]
    print >> sys.stderr, "made word list"
    tree = BKTree(words)
    print >> sys.stderr, "made tree"
    print tree.nodes

    for i in range(6):
        print i, "ddd", [(w.word, w.info ) for w in tree.find("ddd", i)]

def test_speed(dict_file='/usr/share/dict/words', n=1000):
    import re
    print >>sys.stderr, "timing..."
    words = [x.strip() for x in open(dict_file) if not re.search("\W", x.strip()) and len(x.strip()) > 2][:n]
    t = time.time()
    for i,ww in enumerate(words):
        w1 = words[i]
        d = None
        for j,ww in enumerate(words):
            w2 = words[j]
            if w1 == w2: continue
            d = edit_distance(w1, w2, 1000)
    print >> sys.stderr, "\n%i iterations in %.3f seconds of %i-sized tree\n" % (n**2, time.time() - t, n)


def test_tree_dict(dict_file='/usr/share/dict/words'):
    import os
    import re
    if not os.path.exists(dict_file): return False

    words = [x.strip() for x in open(dict_file) if not re.search("\W", x.strip()) and len(x.strip()) > 2]

    t0 = time.time()
    print >>sys.stderr, "\ncreating tree with %i words..." % (len(words),)
    bt = BKTree(words)
    t1 = time.time()
    print >>sys.stderr, "time to create tree:", t1 - t0

    find_words = ("word", "frank", "puddle", "alphabet", "fandango", "puzzle", "spectacular")
    find_dists = (1, 1, 1, 3, 3, 1, 3)

    for fd, fw in zip(find_dists, find_words):
        for i in range(1, 10):
            found = bt.find(fw, i)
            if len(found) > 4:
                print fw, "(" + str(i) + ") :", [w.word for w in found]
                assert fd == i
                break

    t2 = time.time()
    print >>sys.stderr, "time to search tree:", t2 - t1
    #print bt.nodes

    
    for w in words[:100]:
        bt.find(w, 2)
    print >>sys.stderr, "time to search tree for all vs. all at distance == 2:", time.time() - t2

def all():
    test_edit_distance()
    test_tree()
    test_speed()
    test_tree_dict()

if __name__ == "__main__":
    test_same_words()
