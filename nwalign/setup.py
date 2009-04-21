from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Distutils import build_ext

version = '0.1'
import numpy
np_include = numpy.get_include()

setup(name='nwalign',
      version=version,
      description="",
      long_description="""\
""",
      ext_modules=[ Extension("nwalign",
                      sources=["nwalign.pyx"],
                      include_dirs=[np_include, ","])],
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bio',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='BSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points= {
          # -*- Entry points: -*-
          'console_scripts': ['nwalign = nwalign:main']
          }
      )
