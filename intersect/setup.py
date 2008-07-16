#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
#

#print 'remember to cython --convert-range intersection.pyx before  setup.py build_ext --inplace'
import os
os.system('cython --convert-range intersection.pyx')

setup( name = 'intersection'
        , ext_modules=[ Extension("intersection", sources=["intersection.pyx"])]
        , cmdclass = {'build_ext': build_ext}
  )

