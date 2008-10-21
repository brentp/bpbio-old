"""
    >>> from interval_hash import dumps, loads

    >>> loads(dumps(0, 32))
    (0, 32)

    >>> loads(dumps(1025, 2048))
    (1024, 2048)

    >>> loads(dumps(1024, 2048))
    (0, 2048)

it cant save the exact interval:

    >>> loads(dumps(1023, 2049))
    (0, 4096)

    >>> b12 = dumps(1, 2)
    >>> b34 = dumps(3, 4)
    >>> b12, b34
    ('00000000000000000000000', '000000000000000000000011')


    >>> b12 < b34
    True


But in most cases that wont work as the precision gets lost:
    >>> b47 = dumps(4, 7)
    >>> b89 = dumps(8, 9)

    >>> b47 < b89
    False

Must always use <= then compare the intervals directly
    >>> a = dumps(40000, 70000)
    >>> b = dumps(80000, 90000)
    >>> a, b
    ('0000000', '0000000111')


    >>> loads(a), loads(b)
    ((0, 131072), (114688, 131072))

"""
def dumps(imin, imax, rng=2**24):
    """take a start (imin) and stop (imax)
    and convert to a representation that indicates
    position and size of the interval. suitable for
    a btree index
        >>> a = dumps(40000, 70000)
        >>> b = dumps(80000, 90000)
        >>> a, b
        ('0000000', '0000000111')

    """

    ilist = []
    i0 = i1 = '0'
    dist = imax - imin
    while i0 == i1 and rng >= dist:
        rng >>= 1
        i0 = imin <= rng and '0' or '1'
        i1 = imax <= rng and '0' or '1'
        ilist.append(i0)

    return ''.join(ilist[:-1])


def loads(hashed, rng=2**24):
    """given a hashed interval, return the smallest interval
    that contains it:
    >>> a, b = ('0000000', '0000000111')
    >>> loads(a), loads(b)
    ((0, 131072), (114688, 131072))

    """


    imin, imax = 0, rng
    
    for io in hashed:
        rng >>= 1
        if io == '0': imax-= rng
        else: imin += rng 
        #print imin, imax
    return (imin, imax)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
