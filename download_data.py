#!/usr/bin/env python3
"""
Download required datasets for predNMD

This script downloads large annotation files needed for predNMD analysis.

Features:
- Static datasets: gnomAD constraint metrics, phyloP conservation scores
- Dynamic Ensembl datasets: GTF annotations, CDS sequences, reference genomes, VEP cache
- Flexible assembly selection: GRCh37/hg19 or GRCh38/hg38
- Flexible Ensembl release version selection
- Progress bars (requires tqdm)
- MD5 checksum verification (when available)

Basic usage:
    # Download static datasets only
    python download_data.py
    
    # Download Ensembl files for GRCh38 (latest release)
    python download_data.py --datasets ensembl-gtf ensembl-cds ensembl-genome
    
    # Download Ensembl VEP cache for GRCh38
    python download_data.py --datasets ensembl-vep
    
    # Download Ensembl files for GRCh37 with specific release
    python download_data.py --datasets ensembl-all --assembly GRCh37 --ensembl-release 87
    
    # List all available datasets
    python download_data.py --list
    
    # Download specific datasets with custom directory
    python download_data.py --datasets gnomad ensembl-gtf ensembl-vep --data-dir /path/to/data

For more examples, run: python download_data.py --help
"""

import os
import sys
import hashlib
import urllib.request
import argparse
from pathlib import Path


# Ensembl release information
LATEST_ENSEMBL_RELEASE = 112  # Update as needed
GRCH37_LAST_RELEASE = 87  # Last release supporting GRCh37 in main FTP

# Ensembl URL templates
ENSEMBL_BASE = 'http://ftp.ensembl.org/pub'
ENSEMBL_GRCH37_BASE = 'http://ftp.ensembl.org/pub/grch37'


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


def get_ensembl_assembly_name(assembly):
    """Convert assembly identifier to Ensembl format"""
    assembly_map = {
        'GRCh37': 'GRCh37',
        'hg19': 'GRCh37',
        'grch37': 'GRCh37',
        'GRCh38': 'GRCh38',
        'hg38': 'GRCh38',
        'grch38': 'GRCh38'
    }
    return assembly_map.get(assembly, assembly)


def get_ensembl_gtf_info(assembly='GRCh38', release=None):
    """Generate Ensembl GTF dataset information"""
    assembly = get_ensembl_assembly_name(assembly)
    
    if release is None:
        release = GRCH37_LAST_RELEASE if assembly == 'GRCh37' else LATEST_ENSEMBL_RELEASE
    
    if assembly == 'GRCh37':
        # GRCh37 uses special URL structure
        url = f"{ENSEMBL_GRCH37_BASE}/release-{release}/gtf/homo_sapiens/Homo_sapiens.{assembly}.{release}.gtf.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.gtf.gz"
    else:
        url = f"{ENSEMBL_BASE}/release-{release}/gtf/homo_sapiens/Homo_sapiens.{assembly}.{release}.gtf.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.gtf.gz"
    
    return {
        'filename': filename,
        'url': url,
        'description': f'Ensembl {assembly} GTF annotation (release {release})',
        'md5': None,
        'size_mb': 50,  # Approximate
    }


def get_ensembl_cds_info(assembly='GRCh38', release=None):
    """Generate Ensembl CDS FASTA dataset information"""
    assembly = get_ensembl_assembly_name(assembly)
    
    if release is None:
        release = GRCH37_LAST_RELEASE if assembly == 'GRCh37' else LATEST_ENSEMBL_RELEASE
    
    if assembly == 'GRCh37':
        url = f"{ENSEMBL_GRCH37_BASE}/release-{release}/fasta/homo_sapiens/cds/Homo_sapiens.{assembly}.cds.all.fa.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.cds.all.fa.gz"
    else:
        url = f"{ENSEMBL_BASE}/release-{release}/fasta/homo_sapiens/cds/Homo_sapiens.{assembly}.cds.all.fa.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.cds.all.fa.gz"
    
    return {
        'filename': filename,
        'url': url,
        'description': f'Ensembl {assembly} CDS sequences (release {release})',
        'md5': None,
        'size_mb': 30,  # Approximate
    }


def get_ensembl_genome_info(assembly='GRCh38', release=None):
    """Generate Ensembl genome FASTA dataset information"""
    assembly = get_ensembl_assembly_name(assembly)
    
    if release is None:
        release = GRCH37_LAST_RELEASE if assembly == 'GRCh37' else LATEST_ENSEMBL_RELEASE
    
    if assembly == 'GRCh37':
        url = f"{ENSEMBL_GRCH37_BASE}/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.{assembly}.dna.primary_assembly.fa.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.dna.primary_assembly.fa.gz"
    else:
        url = f"{ENSEMBL_BASE}/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.{assembly}.dna.primary_assembly.fa.gz"
        filename = f"Homo_sapiens.{assembly}.{release}.dna.primary_assembly.fa.gz"
    
    return {
        'filename': filename,
        'url': url,
        'description': f'Ensembl {assembly} reference genome (release {release})',
        'md5': None,
        'size_mb': 900,  # ~900 MB compressed
    }


def get_ensembl_vep_info(assembly='GRCh38', release=None):
    """Generate Ensembl-VEP cache dataset information"""
    assembly = get_ensembl_assembly_name(assembly)
    
    if release is None:
        release = GRCH37_LAST_RELEASE if assembly == 'GRCh37' else LATEST_ENSEMBL_RELEASE
    
    if assembly == 'GRCh37':
        # GRCh37 uses special URL structure
        url = f"{ENSEMBL_GRCH37_BASE}/release-{release}/variation/indexed_vep_cache/homo_sapiens_vep_{release}_GRCh37.tar.gz"
        filename = f"homo_sapiens_vep_{release}_GRCh37.tar.gz"
    else:
        url = f"{ENSEMBL_BASE}/release-{release}/variation/indexed_vep_cache/homo_sapiens_vep_{release}_{assembly}.tar.gz"
        filename = f"homo_sapiens_vep_{release}_{assembly}.tar.gz"
    
    return {
        'filename': filename,
        'url': url,
        'description': f'Ensembl-VEP cache for {assembly} (release {release})',
        'md5': None,
        'size_mb': 18000,  # ~18 GB compressed (approximate)
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
        description='Download required datasets for predNMD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all datasets
  python download_data.py
  
  # Download only specific datasets
  python download_data.py --datasets gnomad phylop-hg38
  
  # Download Ensembl files for GRCh38
  python download_data.py --datasets ensembl-gtf ensembl-cds ensembl-genome --assembly GRCh38 --ensembl-release 112
  
  # Download Ensembl VEP cache for GRCh38
  python download_data.py --datasets ensembl-vep --assembly GRCh38 --ensembl-release 112
  
  # Download Ensembl files for GRCh37
  python download_data.py --datasets ensembl-gtf ensembl-cds --assembly GRCh37 --ensembl-release 87
  
  # Download all core Ensembl files (GTF, CDS, genome - excludes VEP cache)
  python download_data.py --datasets ensembl-all --assembly GRCh38
  
  # Download Ensembl files including VEP cache
  python download_data.py --datasets ensembl-all ensembl-vep --assembly GRCh38
  
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
        choices=list(DATASETS.keys()) + ['ensembl-gtf', 'ensembl-cds', 'ensembl-genome', 'ensembl-vep', 'ensembl-all'],
        default=None,
        help='Specific datasets to download (default: static datasets only)'
    )
    
    parser.add_argument(
        '--assembly',
        type=str,
        default='GRCh38',
        choices=['GRCh37', 'GRCh38', 'hg19', 'hg38'],
        help='Genome assembly version for Ensembl files (default: GRCh38)'
    )
    
    parser.add_argument(
        '--ensembl-release',
        type=int,
        default=None,
        help=f'Ensembl release version (default: {LATEST_ENSEMBL_RELEASE} for GRCh38, {GRCH37_LAST_RELEASE} for GRCh37)'
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
    
    # Normalize assembly name
    assembly = get_ensembl_assembly_name(args.assembly)
    
    # Set default release if not specified
    if args.ensembl_release is None:
        ensembl_release = GRCH37_LAST_RELEASE if assembly == 'GRCh37' else LATEST_ENSEMBL_RELEASE
    else:
        ensembl_release = args.ensembl_release
    
    # Build dynamic datasets dictionary
    dynamic_datasets = {}
    if args.datasets and any(d.startswith('ensembl-') for d in args.datasets):
        if 'ensembl-gtf' in args.datasets or 'ensembl-all' in args.datasets:
            dynamic_datasets['ensembl-gtf'] = get_ensembl_gtf_info(assembly, ensembl_release)
        if 'ensembl-cds' in args.datasets or 'ensembl-all' in args.datasets:
            dynamic_datasets['ensembl-cds'] = get_ensembl_cds_info(assembly, ensembl_release)
        if 'ensembl-genome' in args.datasets or 'ensembl-all' in args.datasets:
            dynamic_datasets['ensembl-genome'] = get_ensembl_genome_info(assembly, ensembl_release)
        if 'ensembl-vep' in args.datasets:
            dynamic_datasets['ensembl-vep'] = get_ensembl_vep_info(assembly, ensembl_release)
    
    # Merge static and dynamic datasets
    all_datasets = {**DATASETS, **dynamic_datasets}
    
    # List datasets and exit
    if args.list:
        print("\nAvailable datasets:")
        print("=" * 80)
        print("\nStatic datasets:")
        for name, info in DATASETS.items():
            print(f"\n{name}:")
            print(f"  Description: {info['description']}")
            print(f"  Filename: {info['filename']}")
            print(f"  Size: ~{info['size_mb']} MB")
            print(f"  URL: {info['url']}")
        
        print("\n" + "=" * 80)
        print("\nEnsembl datasets (dynamic - specify with --assembly and --ensembl-release):")
        print("\nensembl-gtf:")
        print(f"  Description: Ensembl GTF annotation file")
        print(f"  Example: Homo_sapiens.GRCh38.112.gtf.gz")
        print(f"  Size: ~50 MB")
        print("\nensembl-cds:")
        print(f"  Description: Ensembl CDS sequences (FASTA)")
        print(f"  Example: Homo_sapiens.GRCh38.112.cds.all.fa.gz")
        print(f"  Size: ~30 MB")
        print("\nensembl-genome:")
        print(f"  Description: Ensembl reference genome (FASTA)")
        print(f"  Example: Homo_sapiens.GRCh38.112.dna.primary_assembly.fa.gz")
        print(f"  Size: ~900 MB")
        print("\nensembl-vep:")
        print(f"  Description: Ensembl-VEP cache (indexed)")
        print(f"  Example: homo_sapiens_vep_112_GRCh38.tar.gz")
        print(f"  Size: ~18 GB")
        print("\nensembl-all:")
        print(f"  Description: All core Ensembl files (GTF + CDS + genome, excludes VEP cache)")
        
        return 0
    
    # Determine which datasets to download
    if args.datasets is None:
        # Default: only static datasets
        datasets_to_download = list(DATASETS.keys())
    else:
        # Expand 'ensembl-all' if present
        datasets_to_download = []
        for d in args.datasets:
            if d == 'ensembl-all':
                datasets_to_download.extend(['ensembl-gtf', 'ensembl-cds', 'ensembl-genome'])
            else:
                datasets_to_download.append(d)
    
    # Filter to only available datasets
    datasets_to_download = [d for d in datasets_to_download if d in all_datasets]
    
    if not datasets_to_download:
        print("Error: No valid datasets specified")
        return 1
    
    # Determine data directory
    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = get_data_directory()
    
    # Create data directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("predNMD Data Download")
    print("=" * 80)
    print(f"\nData directory: {data_dir.absolute()}")
    print(f"Datasets to download: {', '.join(datasets_to_download)}")
    
    # Show Ensembl configuration if applicable
    if any(d.startswith('ensembl-') for d in datasets_to_download):
        print(f"\nEnsembl configuration:")
        print(f"  Assembly: {assembly}")
        print(f"  Release: {ensembl_release}")
    
    # Calculate total size
    total_size = sum(all_datasets[d]['size_mb'] for d in datasets_to_download)
    print(f"\nTotal download size: ~{total_size} MB (~{total_size/1024:.1f} GB)")
    
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
    
    for dataset_name in datasets_to_download:
        dataset_info = all_datasets[dataset_name]
        
        if download_dataset(dataset_name, dataset_info, data_dir, args.force):
            success_count += 1
        else:
            failed_datasets.append(dataset_name)
    
    # Print summary
    print("\n" + "=" * 80)
    print("Download Summary")
    print("=" * 80)
    print(f"Successful: {success_count}/{len(datasets_to_download)}")
    
    if failed_datasets:
        print(f"Failed: {', '.join(failed_datasets)}")
        print("\nSome downloads failed. Please check your internet connection")
        print("and try again. You can also download files manually from:")
        for dataset_name in failed_datasets:
            print(f"  {all_datasets[dataset_name]['url']}")
        return 1
    else:
        print("\n✓ All datasets downloaded successfully!")
        print(f"\nData location: {data_dir.absolute()}")
        
        # Show filenames for Ensembl files
        if any(d.startswith('ensembl-') for d in datasets_to_download):
            print("\nEnsembl files downloaded:")
            for dataset_name in datasets_to_download:
                if dataset_name.startswith('ensembl-'):
                    print(f"  {all_datasets[dataset_name]['filename']}")
        
        print("\nYou can now run predNMD with these annotation files.")
        return 0


if __name__ == '__main__':
    sys.exit(main())
