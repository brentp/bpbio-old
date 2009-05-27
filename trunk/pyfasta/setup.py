from setuptools import setup, find_packages
import sys, os


version = '0.1'

setup(name='pyfasta',
      version=version,
      description="pythonic access to fasta sequence files",
      url="http://code.google.com/p/bpbio/",
      long_description=open('README.txt').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      keywords='bioinformatics blast',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
