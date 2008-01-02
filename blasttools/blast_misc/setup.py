from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
        name="blast_misc"

        ,ext_modules = [
            ("blast_misc", ["blast_misc.pyx"])
         ]
        ,cmd_class = {'build_ext':build_ext}
)


