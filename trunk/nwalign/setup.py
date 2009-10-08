from setuptools import setup, find_packages
from distutils.extension import Extension
#from Cython.Distutils import build_ext

version = '0.1'
import numpy
np_include = numpy.get_include()
import nwalign

setup(name='nwalign',
      version=version,
      description="Needleman-Wunsch global sequence alignment",
      long_description=nwalign.__doc__,
      ext_modules=[ Extension("nwalign",
                              #sources=["nwalign.pyx"],
                      sources=["nwalign.c"],
                      include_dirs=[np_include, "."])],
      keywords='bioinformatics alignment needleman-wunsch',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='BSD',
      packages=find_packages(exclude=[]), # + ['.'],
      test_suite='nose.collector',
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'numpy', 'cython'
      ],
      entry_points= {
          # -*- Entry points: -*-
          'console_scripts': ['nwalign = nwalign:main']
          },
    classifiers   = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering',
        'Topic :: Text Processing'
        ],

)