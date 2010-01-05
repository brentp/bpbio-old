from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='gtpym',
      version=version,
      description="""Extensions to genometools gtpython""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bioinformatics',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      #packages=['gtpym'],
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
          #'gt', 'numpy'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
