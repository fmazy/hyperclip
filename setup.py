#!/usr/bin/env python3
from setuptools import setup

__version__ = "0.1.0"


setup(
    name='hyperclip',
    version=__version__,
    description='Volume of Hypercubes Clipped by Hyperplanes',
    author='François-Rémi Mazy',
    author_email='francois-remi.mazy@inria.fr',
    url='https://gitlab.inria.fr/fmazy/hyperclip',
    packages=[ 'hyperclip', ],
    package_dir={
        'hyperclip' : 'hyperclip',
    },
    install_requires=[
            'numpy>=1.19.2',
        ],

    long_description=open('README.md').read(),

    license="GNU GENERAL PUBLIC LICENSE",
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords='volume hypercube hyperplanes clip',
)
