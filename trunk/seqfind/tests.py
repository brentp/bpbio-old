#
from seqfind import edit_distance
from seqfind import BKTree


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


def test_tree():
    words = ("abcd", "def", "acdf", "gef", "wersdf", "asdfewd"
            , "sdfwe", "sdfwef", "kljd", "oiouoij", "ououoji", "ouonkoj")

    tree = BKTree(words)
    print "made tree"

    for i in range(6):
        print i, "ddd", tree.find("ddd", i)

def test_speed(n=500000):
    print "timing..."
    import time
    t = time.time()
    for i in xrange(n):
        edit_distance('i ehm a gude spehlar', 'i am a good speller', 100)
    print "%i iterations in %f" % (n, time.time() - t)


def test_tree_dict(dict_file='/usr/share/dict/words'):
    import os
    import re
    if not os.path.exists(dict_file): return False
    import time, sys
    import random
    words = [x.strip() for x in open(dict_file) if not re.search("\W", x.strip()) and len(x.strip()) > 2]
    random.shuffle(words)

    t0 = time.time()
    print >>sys.stderr, "creating tree..."
    bt = BKTree(words)
    t1 = time.time()
    print >>sys.stderr, "time to create tree:", t1 - t0
    find_words = ("word", "frank", "puddle", "alphabet", "fandango", "puzzle", "spectacular")

    for fw in find_words:
        for i in range(10):
            found = bt.find(fw, i)
            if len(found) > 4:
                print fw, "(" + str(i) + ") :", found
                break

    t2 = time.time()
    print >>sys.stderr, "time to search tree:", t2 - t1

def all():
    test_edit_distance()
    test_tree()
    test_speed()


