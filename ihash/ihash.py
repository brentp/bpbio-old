"""
    >>> from ihash import dumps, loads

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
    ('0000000000000000000000000', '0000000000000000000000001')


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
    ('000000000', '00000000010')

    >>> loads(a), loads(b)
    ((0, 131072), (65536, 98304))


   >>> loads(dumps(15993504, 15993514, 2**24), 2**24)
   (15993472, 15993536)

   >>> loads(dumps(2, 2))
   (0, 2)
"""
def dumps(imin, imax, rng=2**26):
    """take a start (imin) and stop (imax)
    and convert to a representation that indicates
    position and size of the interval. suitable for
    a btree index
        >>> a = dumps(40000, 70000)
        >>> b = dumps(80000, 90000)
        >>> a, b
        ('000000000', '00000000010')

    """

    ilist = []
    i0 = i1 = '0'
    rng >>= 1
    mid = rng
    while rng > 1:
        rng >>=1
        if imin <= mid and imax <= mid:
            mid -= rng
            ilist.append('0')
        elif imin > mid and imax > mid:
            mid += rng
            ilist.append('1')
        else:
            break
    return ''.join(ilist)

def loads(hashed, rng=2**26):
    """given a hashed interval, return the smallest interval
    that contains it:
    >>> a, b = ('000000000', '000000000111')
    >>> loads(a), loads(b)
    ((0, 131072), (114688, 131072))

    """


    imin, imax = 0, rng
    for io in hashed:
        rng >>= 1
        if io == '0': imax-= rng;
        else: imin += rng; 
    return (imin, imax)

def precision(length, rng=2**26):
    return rng >> length

try:
    from _ihash import dumps, loads
except:
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()

