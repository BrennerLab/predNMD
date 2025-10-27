#!/usr/bin/env python3
"""
Download required datasets for deepNMD

This script downloads large annotation files needed for deepNMD analysis.
Run this after cloning the repository:
    python download_data.py
"""

import os
import sys
import hashlib
import urllib.request
import argparse
from pathlib import Path


# Define dataset URLs and metadata
DATASETS = {
    'gnomad': {
        'filename': 'gnomad.v4.1.constraint_metrics.tsv',
        'url': 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv',
        'description': 'gnomAD v4.1 constraint metrics',
        'md5': None,  # Add checksum if available
        'size_mb': 90,  # Approximate size
    },
    'phylop-hg19': {
        'filename': 'hg19.100way.phyloP100way.bw',
        'url': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw',
        'description': 'phyloP conservation scores for hg19 (BigWig)',
        'md5': None,
        'size_mb': 9500,  # ~9.5GB
    },
    'phylop-hg38': {
        'filename':'hg38.phyloP100way.bw',
        'url':'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw',
        'description':'phyloP conservation scores for hg38 (BigWig)',
        'md5': None,
        'size_mb': 9500,  # ~9.5GB
    }
}


try:
    from tqdm import tqdm
    
    class DownloadProgressBar(tqdm):
        """Progress bar for downloads"""
        def update_to(self, b=1, bsize=1, tsize=None):
            if tsize is not None:
                self.total = tsize
            self.update(b * bsize - self.n)
    
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 checksum of a file"""
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)
    return md5.hexdigest()


def verify_checksum(filepath, expected_md5):
    """Verify file checksum"""
    if expected_md5 is None:
        return True
    
    print(f"  Verifying checksum...")
    actual_md5 = calculate_md5(filepath)
    
    if actual_md5 == expected_md5:
        print(f"  ✓ Checksum verified")
        return True
    else:
        print(f"  ✗ Checksum mismatch!")
        print(f"    Expected: {expected_md5}")
        print(f"    Got: {actual_md5}")
        return False


def download_file(url, output_path, description=None):
    """Download a file with progress bar"""
    print(f"\nDownloading: {description or output_path}")
    print(f"  URL: {url}")
    print(f"  Destination: {output_path}")
    
    try:
        if HAS_TQDM:
            with DownloadProgressBar(unit='B', unit_scale=True, miniters=1,
                                    desc=str(output_path.name)) as t:
                urllib.request.urlretrieve(
                    url, 
                    filename=output_path,
                    reporthook=t.update_to
                )
        else:
            # Simple download without progress bar
            print(f"  Downloading... (install tqdm for progress bar)")
            urllib.request.urlretrieve(url, filename=output_path)
        
        print(f"  ✓ Download complete")
        return True
    except Exception as e:
        print(f"  ✗ Download failed: {e}")
        if output_path.exists():
            output_path.unlink()
        return False


def download_dataset(dataset_name, dataset_info, data_dir, force=False):
    """Download a single dataset"""
    output_path = data_dir / dataset_info['filename']
    
    # Check if file already exists
    if output_path.exists() and not force:
        print(f"\n✓ {dataset_info['filename']} already exists")
        if dataset_info['md5']:
            if verify_checksum(output_path, dataset_info['md5']):
                return True
            else:
                print(f"  Re-downloading due to checksum mismatch...")
        else:
            return True
    
    # Download the file
    if download_file(dataset_info['url'], output_path, dataset_info['description']):
        # Verify checksum if provided
        if dataset_info['md5']:
            if not verify_checksum(output_path, dataset_info['md5']):
                print(f"  Warning: Checksum verification failed")
                return False
        return True
    
    return False


def get_data_directory():
    """Get the data directory path"""
    # Try to find the package directory
    script_dir = Path(__file__).parent.absolute()
    
    # Check if we're in the repo root
    if (script_dir / 'setup.py').exists() or (script_dir / 'pyproject.toml').exists():
        data_dir = script_dir / 'data'
    # Check if we're in a scripts directory
    elif (script_dir.parent / 'setup.py').exists():
        data_dir = script_dir.parent / 'data'
    else:
        # Default to current directory
        data_dir = script_dir / 'data'
    
    return data_dir


def main():
    parser = argparse.ArgumentParser(
        description='Download required datasets for deepNMD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all datasets
  python download_data.py
  
  # Download only specific datasets
  python download_data.py --datasets gnomad phylop
  
  # Force re-download even if files exist
  python download_data.py --force
  
  # Specify custom data directory
  python download_data.py --data-dir /path/to/data
        """
    )
    
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=None,
        help='Directory to store downloaded data (default: ./data)'
    )
    
    parser.add_argument(
        '--datasets',
        nargs='+',
        choices=list(DATASETS.keys()),
        default=list(DATASETS.keys()),
        help='Specific datasets to download (default: all)'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download even if files exist'
    )
    
    parser.add_argument(
        '--list',
        action='store_true',
        help='List available datasets and exit'
    )
    
    args = parser.parse_args()
    
    # List datasets and exit
    if args.list:
        print("\nAvailable datasets:")
        print("=" * 80)
        for name, info in DATASETS.items():
            print(f"\n{name}:")
            print(f"  Description: {info['description']}")
            print(f"  Filename: {info['filename']}")
            print(f"  Size: ~{info['size_mb']} MB")
            print(f"  URL: {info['url']}")
        return 0
    
    # Determine data directory
    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = get_data_directory()
    
    # Create data directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("deepNMD Data Download")
    print("=" * 80)
    print(f"\nData directory: {data_dir.absolute()}")
    print(f"Datasets to download: {', '.join(args.datasets)}")
    
    # Calculate total size
    total_size = sum(DATASETS[d]['size_mb'] for d in args.datasets)
    print(f"Total download size: ~{total_size} MB (~{total_size/1024:.1f} GB)")
    
    # Confirm with user
    if not args.force:
        response = input("\nProceed with download? [y/N]: ")
        if response.lower() != 'y':
            print("Download cancelled.")
            return 1
    
    # Download datasets
    print("\nStarting downloads...")
    print("-" * 80)
    
    success_count = 0
    failed_datasets = []
    
    for dataset_name in args.datasets:
        dataset_info = DATASETS[dataset_name]
        
        if download_dataset(dataset_name, dataset_info, data_dir, args.force):
            success_count += 1
        else:
            failed_datasets.append(dataset_name)
    
    # Print summary
    print("\n" + "=" * 80)
    print("Download Summary")
    print("=" * 80)
    print(f"Successful: {success_count}/{len(args.datasets)}")
    
    if failed_datasets:
        print(f"Failed: {', '.join(failed_datasets)}")
        print("\nSome downloads failed. Please check your internet connection")
        print("and try again. You can also download files manually from:")
        for dataset_name in failed_datasets:
            print(f"  {DATASETS[dataset_name]['url']}")
        return 1
    else:
        print("\n✓ All datasets downloaded successfully!")
        print(f"\nData location: {data_dir.absolute()}")
        print("\nYou can now run deepNMD with these annotation files.")
        return 0


if __name__ == '__main__':
    sys.exit(main())
