from setuptools import setup


setup(name          = 'ncsv',
      version       = '0.1.0',
      description   = 'easier csv reading/writing from numpy',
      license       = 'MIT',
      keywords      = 'numpy csv',
      author        = 'Brent Pedersen',
      author_email  = 'bpederse@gmail.com',
      url   = '',
      packages      = [''],
      requires = ['numpy'],
      test_suite = 'tests',
      classifiers   = [
        'Development Status :: 4 - Beta'
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
        ],
)
