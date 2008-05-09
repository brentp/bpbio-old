
#from distutils.core import setup 
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
#

print 'remember to cython seqfind.pyx  setup.py install'

setup( name = 'seqfind'
        , ext_modules=[ Extension("seqfind", sources=["seqfind.pyx"] 
           # , language="c++"
           )]
        , cmdclass = {'build_ext': build_ext}
        , test_suite = 'nose.collector'
  )

