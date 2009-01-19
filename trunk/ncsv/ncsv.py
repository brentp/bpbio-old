import numpy as np

__all__ = ("savetxt", "loadtxt")

def savetxt(A, filename, sep="\t", names=None):
    fields = _tuple_dtype(A)
    if names is not None:
        fields = [(name, dt) for name, (old_name, dt) in zip(names, fields)]
    # dtype.str is e.g. '<i4' so get just 'i4'
    header = sep.join(["%s||%s" % (fname, dt.str[1:]) for fname, dt in fields])

    if hasattr(filename, 'seek'):
        fh = filename
        fh.seek(0)
    else:
        fh = open(filename, 'wb')
    fh.write('#' + header + '\n')
    fmt = _fmt_str(A)
    np.savetxt(fh, A, fmt, sep)
    

def _fmt_str(A):
    """convert the numpy dtypes to a python format string
     with %s, %f, %i, etc""" 
    types = [td[1].str[1].lower() for td in _tuple_dtype(A)]
    return ["%" + f for f in types]


def _tuple_dtype(A):
    dt = A.dtype.fields
    t = [(field_name, dtype) for field_name, (dtype, s) in dt.iteritems()]
    return t

def _guess_dtype(filename, sep, n=1000):
    """given a csv filename
    the separator : "\t", "," ... etc.
    and the number of lines to check
    return a dtype suitable for use in the numpy array
    constructor
    """

    fh = open(filename, 'rb')
    l = fh.readline().rstrip()
    names = []
    while l.startswith('#') or not l:
        if l.startswith('#') and names == []:
            names = l.split(sep)
        l = fh.readline().rstrip()

    cols = l.split(sep)
    ncols = len(cols)
    guess = []
    for col in cols:
        guess.append(_get_base_type(col))


    assert len(guess) == len(cols)

def _get_base_type(col):
    # todo. determine the width for float, int
    # 32bit in max is 2147483647
    try:
        int(col)
        return 'i8'
    except:
        try:
            float(col)
            return 'f8'
        except:
            str(col)
            return 'S' + len(col)
    raise

def loadtxt(filename, dtype=None, skiprows=0, sep="\t"):
    if hasattr(filename, 'seek'):
        fh = filename
        fh.seek(0)
    else:
        fh = open(filename, 'rb')

    if dtype is None:
        header = fh.readline().rstrip('#')
        # it has a header.
        if header[0] == '#':
            names, types = zip(*[x.split("||") for x in header.split(sep)])
            dtype = dict(names=names, formats=types)
        else:
            # no header. guess dtype.
            dtype = _guess_dtype(filename, sep)

    return np.loadtxt(fh, dtype=dtype, skiprows=skiprows, delimiter=sep)
