#!/usr/bin/env python
from setuptools import setup
setup(
    name = "cyclopy",
    version = "0.0.1",
    packages = ['cyclopy'],
    install_requires = ['numpy','scipy','mtspec'],

    # metadata for upload to PyPI
    author = "Luca Penasa",
    author_email = "luca.penasa@gmail.com",
    description = "Cyclostratigraphic analysis made easy in python",
    #license = "BSD",
    #keywords = "",
    #url = "",   
    classifiers =
        [ 'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: C',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Operating System :: OS Independent',
          ],
)
