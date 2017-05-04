#!/usr/bin/env python

from setuptools import setup
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
# import sys



setup(
    name = "cyclopy",
    version = "0.0.1",
    packages = ['cyclopy', 'cyclopy.Elements', 'cyclopy.External', 'cyclopy.GeologyHelpers','cyclopy.Data', 'cyclopy.NoiseModeling',
                'cyclopy.NumericalMethods', 'cyclopy.Orbitals', 'cyclopy.Utils', 'cyclopy.tests'],
    install_requires = ['numpy','scipy','mtspec','matplotlib'],

    # metadata for upload to PyPI
    author = "Luca Penasa",
    author_email = "luca.penasa@gmail.com",
    description = "Cyclostratigraphic analysis ported to python",
#    license = "GPL",
    keywords = ["cyclostratigraphy", "climate", "time-series", "MTM", "spectral analysis", "geology"],
    url = "https://github.com/luca-penasa/cyclopy",   
    classifiers =
        [ 'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
 #         'License :: OSI Approved :: BSD License',
          'Programming Language :: C',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Operating System :: OS Independent',
          ],
)
