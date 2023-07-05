#!/usr/bin/env python

from setuptools import setup, find_packages

install_requires = [
    'biopython>=1.7.4',
    'gffutils>=0.9',
    'json-source-map>=1.0.1',
    'jsonschema>=3.2.0',
    'pandas>=2.0.0',
    'flask>=2.0.0'
]

setup(
    name='midoe',
    version='1.0',
    packages=find_packages(),
    install_requires=install_requires,
    python_requires='>=3.9',
    classifiers=[
        "Programming Language :: Python :: 3 :: Only"
    ]
)
