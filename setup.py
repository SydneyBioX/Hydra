#!/usr/bin/python

##############################################

# Manoj M Wagle (USydney, MIT CSAIL)

##############################################

import subprocess
import sys
from setuptools import setup, find_namespace_packages

setup(
    name="hydra-tools",
    version="0.1.0",
    packages=find_namespace_packages(),
    install_requires=[
        'anndata',
        'captum',
        'h5py',
        'libpng-bins',
        'torch',
        'torchvision',
        'torchaudio',
        'pandas',
        'scipy',
        'scikit-learn',
        'scanpy',
        'tqdm',
        'numpy',
        'numba'
    ],
    entry_points={
        'console_scripts': [
            'hydra=hydra.Hydra:main',
        ],
    },
    include_package_data=True,
    description='Interpretable deep generative tool for single-cell omics',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/SydneyBioX/Hydra',
    author='Manoj M Wagle',
    author_email='mwag8019@uni.sydney.edu.au',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8'
    ],
)
