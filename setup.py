#!/usr/bin/env python
# coding=utf-8

from distutils.core import setup

# noinspection PyArgumentList
setup(
    name='fatgraph',
    version='1.0.1',
    author='Yuki Koyanagi',
    author_email='yuki@math.au.dk',
    packages=['fatgraph'],
    #url='https://github.com/eseraygun/python-alignment',
    license='BSD 3-Clause License',
    requires=['permutation',],
    package_dir={'fatgraph': 'fatgraph'},
    scripts=['scripts/compute.py',]
)
