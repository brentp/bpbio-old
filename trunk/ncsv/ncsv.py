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


def loadtxt(filename, dtype=None, skiprows=0, sep="\t"):
    if hasattr(filename, 'seek'):
        fh = filename
        fh.seek(0)
    else:
        fh = open(filename, 'rb')

    if dtype is None:
        header = fh.readline().strip('#\n')
        names, types = zip(*[x.split("||") for x in header.split(sep)])
        dtype = dict(names=names, formats=types)

    return np.loadtxt(fh, dtype=dtype, skiprows=skiprows, delimiter=sep)
