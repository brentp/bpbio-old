from setuptools import setup

import simpletable
version = '0.1'

setup(name          = 'simpletable',
      version       = version,
      description   = 'wrapper around pytables/hd5f to simplify using structured data',
      license       = 'BSD',
      keywords      = 'hdf5 pytables tables numpy',
      author        = 'Brent Pedersen',
      author_email  = 'bpederse@gmail.com',
      url   = 'http://bpbio.googlecode.com/',
      download_url = 'http://bpbio.googlecode.com/svn/trunk/simpletable/',
      long_description = simpletable.__doc__,
      install_requires = ['tables'],
      zip_safe=False,
      py_modules = ['simpletable'],
      classifiers   = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Database'
        ],
)
