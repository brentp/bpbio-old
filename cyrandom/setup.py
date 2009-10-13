from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
  name = 'cyrandom',
  ext_modules=[
    Extension("cyrandom",
              sources=["cyrandom.pyx"],
              include_dirs=["."]),
    ],
  cmdclass = {'build_ext': build_ext},

)
