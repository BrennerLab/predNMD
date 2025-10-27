"""
Setup configuration for deepNMD with data download support
"""

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
import subprocess
import sys
import os


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)
        
        # Get the directory where setup.py is located
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        download_script = os.path.join(setup_dir, 'download_data.py')
        
        # Check if user wants to download data
        if os.environ.get('DEEPNMD_DOWNLOAD_DATA', '').lower() in ['true', '1', 'yes']:
            print("\n" + "="*80)
            print("Running post-install data download...")
            print(f"Setup directory: {setup_dir}")
            print(f"Download script: {download_script}")
            print("="*80)
            
            if not os.path.exists(download_script):
                print(f"\nWarning: download_data.py not found at {download_script}")
                print("Skipping data download.")
            else:
                try:
                    # Run the download script with explicit working directory
                    subprocess.check_call(
                        [sys.executable, download_script],
                        cwd=setup_dir,
                        env=os.environ.copy()
                    )
                    print("\nData download completed successfully!")
                except subprocess.CalledProcessError as e:
                    print(f"\nWarning: Data download failed with error: {e}")
                    print(f"You can run 'python {download_script}' manually later.")
                except Exception as e:
                    print(f"\nWarning: Unexpected error during data download: {e}")
                    print(f"You can run 'python {download_script}' manually later.")
        else:
            print("\n" + "="*80)
            print("deepNMD installation complete!")
            print("="*80)
            print("\nTo download required reference data, run:")
            print(f"  python {download_script if os.path.exists(download_script) else 'download_data.py'}")
            print("\nOr set DEEPNMD_DOWNLOAD_DATA=1 before installing:")
            print("  DEEPNMD_DOWNLOAD_DATA=1 pip install -e .")
            print("="*80)


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        
        # Get the directory where setup.py is located
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        download_script = os.path.join(setup_dir, 'download_data.py')
        
        # Check if user wants to download data
        if os.environ.get('DEEPNMD_DOWNLOAD_DATA', '').lower() in ['true', '1', 'yes']:
            print("\n" + "="*80)
            print("Running post-install data download...")
            print(f"Setup directory: {setup_dir}")
            print(f"Download script: {download_script}")
            print("="*80)
            
            if not os.path.exists(download_script):
                print(f"\nWarning: download_data.py not found at {download_script}")
                print("Skipping data download.")
            else:
                try:
                    # Run the download script with explicit working directory
                    subprocess.check_call(
                        [sys.executable, download_script],
                        cwd=setup_dir,
                        env=os.environ.copy()
                    )
                    print("\nData download completed successfully!")
                except subprocess.CalledProcessError as e:
                    print(f"\nWarning: Data download failed with error: {e}")
                    print(f"You can run 'python {download_script}' manually later.")
                except Exception as e:
                    print(f"\nWarning: Unexpected error during data download: {e}")
                    print(f"You can run 'python {download_script}' manually later.")
        else:
            print("\n" + "="*80)
            print("deepNMD installation complete!")
            print("="*80)
            print("\nTo download required reference data, run:")
            print(f"  python {download_script if os.path.exists(download_script) else 'download_data.py'}")
            print("\nOr reinstall with:")
            print("  DEEPNMD_DOWNLOAD_DATA=1 pip install -e .")
            print("="*80)


# Read requirements
with open('requirements.txt') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read long description
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='deepnmd',
    version='1.0.0',
    description='Deep learning-based prediction of NMD escape mechanisms',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Yaqi Su',
    author_email='yaqisu@berkeley.edu',
    url='https://github.com/yaqisu/NMDPred', #TODO: update URL if needed
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
