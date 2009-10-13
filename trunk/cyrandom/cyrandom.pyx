
cdef extern from "stdlib.h":
    void c_srand "srand" (unsigned int)
    int c_rand "rand" ()
    int RAND_MAX


def srand(unsigned int seed):
    c_srand(seed)

seed = srand


cdef inline double rand_float():
    return c_rand() / (RAND_MAX + 0.0)


def random():
    return rand_float()

cpdef int randrange(int min, int max):
    return randint(min, max - 1)

cpdef int randint(int min, int max):
    return <int>(min + ((1.0 + max - min) * rand_float()))


def test():
    import collections
    d = collections.defaultdict(int)
    cdef int i
    for i in range(40000):
        d[randint(1, 4)] += 1
        assert 0.0 <= rand_float() <= 1.0

    assert sorted(d.keys()) == [1, 2, 3, 4]

    # uniformly distributed so should
    # expect around 10,000, should be between
    # 9600 and 10400 "most" of the time.
    # >>> from scipy.stats import binom_test
    # >>> binom_test(10400, 40000, 0.25)
    # >>> 4.3076437759654466e-06

    for v in d.values():
        assert 9600 < v < 10400, v

    # after seeding with same value, should get same result.
    c_srand(1)
    r1 = random()
    c_srand(1)
    r2 = random()
    assert r1 == r2, (r1, r2)
    # next val should be different.
    assert r2 != random()

    # different seed should have different results.
    c_srand(2)
    assert random() != r1
