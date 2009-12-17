from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

version = '0.0'

setup(name='cdagline',
      version=version,
      description="",
      long_description="""\
""",
      cmdclass= {'build_ext': build_ext },
      ext_modules=[ Extension("cdagline",
                      sources=["cdagline.pyx"],)],
      keywords='bio',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='BSD',
      test_suite='nose.collector',
      zip_safe=False,
  )
