from setuptools import setup, find_packages
import sys, os


version = '0.1'

setup(name='pyfasta',
      version=version,
      description="a convenience class for dealing with fasta files in python"
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bioinformatics blast',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
