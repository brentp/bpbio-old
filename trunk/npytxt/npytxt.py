import numpy as np
import os, sys



def loadtxt(fname, npy=True, **kwargs):
    assert os.path.exists(fname)

    npyfile = fname + ".npy"
    if npy and os.path.exists(npyfile) and \
             os.stat(npyfile).st_mtime >= os.stat(fname).st_mtime:
        return np.load(npyfile)

    if not "dtype" in kwargs:  kwargs['dtype'] = None
    if not "names" in kwargs:  kwargs['names'] = True
    A = np.genfromtxt(fname, **kwargs)
    np.save(npyfile, A)
    return A

loadtxt.__doc__ = """\
    a modified version of numpy.loadtxt which will also save a
    pickled version of the file which will be used for faster
    loading after the first time.
    """ + np.loadtxt.__doc__


def savetxt(fname, A, delimiter="\t", names=None):
    npyfile = fname + ".npy"
    names = names or A.dtype.names

    if hasattr(fname, 'seek'):
        fh = fname
        fh.seek(0)
    else:
        fh = open(fname, 'wb')

    if names:
        header = "#" + delimiter.join(names)
        fh.write(header + '\n')
        types = [x[0].str for x in A.dtype.fields.values()]
        # >i4 becomes '%i'
        fmt = ['%' + dt[1].lower() for dt in types]
    else:
        fmt = A.dtype.str[1].lower()

    np.savetxt(fh, A, fmt, delimiter)

    old_names = None
    if names == A.dtype.names:
        np.save(npyfile, A)
        return

    old_names = A.dtype.names
    A.dtype.names = names
    np.save(npyfile, A)
    A.dtype.names = old_names

savetxt.__doc__ = """\
    an extension to numpy.savetxt which will attempt to write a header row
    based on the dtype of the array if a list of names is not given.
    It will save a pickled file which (in combination with loadtxt in 
    this module, will speed loading
    """ + np.savetxt.__doc__

