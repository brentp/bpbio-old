import numpy as np
import os



def loadtxt(fname, **kwargs):
    pklfile = fname + ".nkl"
    if os.path.exists(pklfile) and \
           os.stat(pklfile).st_mtime >= os.stat(fname).st_mtime:
         return np.load(pklfile)
   
    A = np.genfromtxt(fname, **kwargs)
    A.dump(pklfile)
    return A

loadtxt.__doc__ = """\
    a modified version of numpy.loadtxt which will also save a
    pickled version of the file which will be used for faster
    loading after the first time.
    """ + np.loadtxt.__doc__


def savetxt(A, fname, delimiter="\t", names=None):
    pklfile = fname + ".nkl"
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
    # it has a dictionary interface, but seems to preserve order


    np.savetxt(fh, A, fmt, delimiter)
    old_names = None
    if names == A.dtype.names:
        A.dump(pklfile)
        return

    old_names = A.dtype.names
    A.dtype.names = names
    A.dump(pklfile)
    A.dtype.names = old_names

savetxt.__doc__ = """\
    an extension to numpy.savetxt which will attempt to write a header row
    based on the dtype of the array if a list of names is not given.
    It will save a pickled file which (in combination with loadtxt in 
    this module, will speed loading
    """ + np.savetxt.__doc__



