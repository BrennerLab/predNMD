"""
Setup configuration for deepNMD with data download support
"""

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
import subprocess
import sys
import os


# Read requirements
with open('requirements.txt') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read long description
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='nmd',
    version='1.0.0',
    description='Machine learning-based prediction of NMD escape mechanisms',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Yaqi Su',
    author_email='yaqisu@berkeley.edu',
    url='https://github.com/yaqisu/NMD', #TODO: update URL if needed
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'deepnmd': [
            'data/*',
            'models/*.joblib',
            'models/*.json',
        ],
    },
    install_requires=requirements,
    extras_require={
        'dev': [
            'pytest>=7.0',
            'pytest-cov>=4.0',
            'black>=23.0',
            'flake8>=6.0',
            'mypy>=1.0',
        ],
        'download': [
            'tqdm>=4.65',  # For download progress bars
        ],
    },
    entry_points={
        'console_scripts': [
            'deepnmd=deepnmd.cli:main',
            'deepnmd-download-data=download_data:main',
        ],
    },
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='genomics NMD bioinformatics machine-learning',
)
