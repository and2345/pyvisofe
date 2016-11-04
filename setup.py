#!/usr/bin/python2

import codecs
import os
import re
import sys

from setuptools import setup, find_packages

# define some helper functions
# ----------------------------
here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    """
    Build an absolute path from parts and return the contents of the
    resulting file assuming utf-8 encoding.
    """

    with codecs.open(os.path.join(here, *parts), 'rb', 'utf-8') as f:
        return f.read()

meta_path = os.path.join("pyvisofe", "__init__.py")
meta_file = read(meta_path)

def find_meta(meta):
    """
    Extract __*meta*__ from the meta data file.
    """

    meta_match = re.search(pattern=r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
                           string=meta_file,
                           flags=re.M)

    if meta_match:
        return meta_match.group(1)

    raise RuntimeError("Unable to find __{meta}__ string!".format(meta=meta))

# set meta data
# -------------
meta_data = dict(author = "Andreas Kunze",
                 author_email = "andreas.kunze@mailbox.tu-dresden.de",
                 classifiers = ['Environment :: Console',
                                'Development Status :: 3 - Alpha',
                                'Intended Audience :: Science/Research',
                                'Intended Audience :: Developers',
                                'License :: OSI Approved :: BSD License',
                                'Natural Language :: English',
                                'Operating System :: OS Independent',
                                'Programming Language :: Python',
                                'Topic :: Scientific/Engineering',
                                'Topic :: Software Development :: Libraries :: Python Modules'],
                 description = "Scientific visualization software",
                 install_requires = ['numpy>=1.10.4',
                                     'scipy>=0.16.1',
                                     'vispy>=0.5.0.dev0',
                                     'ipython'],
                 license = "BSD",
                 long_description = read("README.rst"),
                 name = "pyvisofe",
                 packages = find_packages(),
                 platforms = 'any',
                 url = "https://github.com/and2345/pyvisofe",
                 version = find_meta('version'),
                 zip_safe=False)

if __name__ == "__main__":
    setup(**meta_data)
