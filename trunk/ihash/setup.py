from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(name = 'ihash'
        , packages=['ihash']
        , ext_modules=[ Extension("_ihash", sources=["ihash/_ihash.pyx"])]
        , cmdclass = {'build_ext': build_ext}
  )

