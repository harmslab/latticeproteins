try:
    from setuptools import setup
except:
    from distutils.core import setup

from distutils.extension import Extension

contactlooper = Extension('latticeproteins.contactlooper', sources = ['latticeproteins/contactlooper.c'])

setup (name = 'latticeproteins',
       fullname = 'Lattice Protein Simulation Package',
       version = '0.1',
       author = 'Jesse D. Bloom',
       author_email = 'jbloom@fhcrc.rog',
       description = 'Code for lattice protein simulations.',
       platforms = 'Tested on Mac OS X and Linux',
       packages = ['latticeproteins'],
       ext_modules = [contactlooper]
)
