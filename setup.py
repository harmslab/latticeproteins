try:
    from setuptools import setup
except:
    from distutils.core import setup

from distutils.extension import Extension

contactlooper = Extension('latticeproteins.contactlooper', sources = ['latticeproteins/contactlooper.c'])

setup(
    name = 'latticeproteins',
    fullname = 'Lattice Protein Simulation Package',
    version = '0.2',
    author = 'Jesse D. Bloom and Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    description = 'Code for lattice protein simulations.',
    url='https://github.com/zsailer/latticeproteins',
    packages = ['latticeproteins'],
    ext_modules = [contactlooper],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
