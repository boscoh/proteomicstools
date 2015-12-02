#!/usr/bin/env python
from setuptools import setup

description = \
"""
A Python module to interface with the proteomics analysis
pipeline: the TPP.
"""

setup(
    name='tppit',
    version='0.9',
    author='Bosco Ho',
    author_email='boscoh@gmail.com',
    url='http://github.com/boscoh/tppit',
    description='module to interface with proteomics TPP package',
    long_description=description,
    license='GPLv3',
    install_requires=['jinja2', 'envoy'],
    py_modules=['tppit',],
    scripts=[],
)