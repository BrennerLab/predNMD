#!/usr/bin/env python3

"""
VEP Annotation and Filtering Pipeline

This script:
1. Annotates VCF with VEP
2. Filters for stop_gained and frameshift variants
3. Outputs filtered VCF file
"""

import argparse
import subprocess
import sys
import os
import gzip
from pathlib import Path


def run_command(cmd, description="Running command"):
    """
    Run a shell command and handle errors.
    
    Args:
        cmd: Command to run (list or string)
        description: Description for logging
    """
    print(f"{description}...")
    print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    try:
        if isinstance(cmd, list):
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        else:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        
        print(f"✓ {description} completed successfully")
        return result
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Error in {description}")
        print(f"Return code: {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise


def run_vep_annotation(vep_path, cache_dir, input_vcf, output_prefix, assembly='GRCh38', threads=32):
    """
    Run VEP annotation on input VCF.
    
    Args:
        vep_path: Path to VEP executable
        cache_dir: VEP cache directory
        input_vcf: Input VCF file (can be .vcf or .vcf.gz)
        output_prefix: Output prefix for files
        assembly: Genome assembly 
        threads: Number of threads for VEP
        
    Returns:
        Path to VEP output file
    """
    vep_output = f"{output_prefix}.vep.vcf"
    
    cmd = [
        vep_path,
        '--cache',
        '--dir_cache', cache_dir,
        '-i', input_vcf,
        '-o', vep_output,
        '--offline',
        '--species', 'homo_sapiens',
        '--assembly', assembly,
        '--vcf',
        '--af_gnomad',
        '--canonical',
        '--biotype',
        '--force_overwrite',
        '--fork', str(threads)
    ]
    
    run_command(cmd, "VEP annotation")
    
    if not os.path.exists(vep_output):
        raise FileNotFoundError(f"VEP output file not found: {vep_output}")
    
    return vep_output


def filter_vep_output(vep_file, output_file):
    """
    Filter VEP output (VCF format) for stop_gained and frameshift variants.
    
    Args:
        vep_file: Path to VEP output file (VCF format)
        output_file: Path to filtered output file
        
    Returns:
        Number of variants found and written
    """
    print("Filtering VEP output for stop_gained and frameshift variants...")
    
    # Target consequences to keep
    target_consequences = {'stop_gained', 'frameshift_variant'}
    
    variants_written = 0
    csq_format = None
    consequence_counts = {}
    
    with open(vep_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Write all header lines and extract CSQ format
            if line.startswith('#'):
                outfile.write(line)
                
                # Extract CSQ format from header
                if line.startswith('##INFO=<ID=CSQ'):
                    # Extract format from Description
                    # Format: Allele|Consequence|IMPACT|...
                    format_start = line.find('Format: ')
                    if format_start != -1:
                        format_str = line[format_start + 8:]
                        format_end = format_str.find('"')
                        if format_end != -1:
                            format_str = format_str[:format_end]
                            csq_format = format_str.split('|')
                continue
            
            # Skip if no CSQ format found
            if not csq_format:
                continue
            
            # Parse variant line
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            info = fields[7]
            
            # Extract CSQ from INFO field
            csq_data = None
            for info_field in info.split(';'):
                if info_field.startswith('CSQ='):
                    csq_data = info_field[4:]  # Remove 'CSQ='
                    break
            
            if not csq_data:
                continue
            
            # Check if any annotation has target consequences
            has_target_consequence = False
            
            # Split by comma to get individual annotations (transcripts)
            annotations = csq_data.split(',')
            
            for annotation in annotations:
                # Split by pipe to get fields
                ann_fields = annotation.split('|')
                
                if len(ann_fields) != len(csq_format):
                    continue
                
                # Create annotation dictionary
                ann_dict = dict(zip(csq_format, ann_fields))
                
                # Check consequence field
                consequence = ann_dict.get('Consequence', '')
                
                for target in target_consequences:
                    if target in consequence:
                        has_target_consequence = True
                        consequence_counts[target] = consequence_counts.get(target, 0) + 1
            
            # Write variant if it has target consequence
            if has_target_consequence:
                outfile.write(line)
                variants_written += 1
    
    print(f"Found and wrote {variants_written} variants with target consequences")
    
    if consequence_counts:
        print("Breakdown by consequence type (transcript-level annotations):")
        for consequence, count in sorted(consequence_counts.items()):
            print(f"  {consequence}: {count}")
    
    return variants_written


def main():
    parser = argparse.ArgumentParser(
        description="VEP annotation and filtering pipeline for stop_gained and frameshift variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with GRCh38 (default)
  python vepAnnot_filter.py input.vcf output_prefix \\
    --vep-path ./vep --cache-dir /path/to/.vep
  
  # With GRCh37 assembly
  python vepAnnot_filter.py input.vcf.gz output_prefix \\
    --vep-path ./vep --cache-dir /path/to/.vep --assembly GRCh37
        """
    )
    
    # Required arguments
    parser.add_argument('input_vcf', help='Input VCF file')
    parser.add_argument('output_prefix', help='Output prefix for generated files')
    
    # VEP arguments
    parser.add_argument('--vep-path', required=True, help='Path to VEP executable')
    parser.add_argument('--cache-dir', required=True, help='VEP cache directory')
    parser.add_argument('--assembly', default='GRCh38', help='Genome assembly (default: GRCh38)')
    parser.add_argument('--threads', type=int, default=32, help='Number of threads for VEP (default: 32)')
    
    # Optional arguments
    parser.add_argument('--keep-full-vep', action='store_true', 
                       help='Keep the full VEP output file (before filtering)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.input_vcf):
        print(f"Error: Input VCF not found: {args.input_vcf}")
        sys.exit(1)
    
    # Check if input is gzipped and readable
    try:
        if args.input_vcf.endswith('.gz'):
            with gzip.open(args.input_vcf, 'rt') as f:
                next(f)  # Try to read first line
            print(f"Input VCF is gzipped: {args.input_vcf}")
        else:
            with open(args.input_vcf, 'r') as f:
                next(f)  # Try to read first line
            print(f"Input VCF is uncompressed: {args.input_vcf}")
    except Exception as e:
        print(f"Error: Cannot read input VCF file: {e}")
        sys.exit(1)
    
    # Check if VEP executable exists
    if not os.path.exists(args.vep_path):
        print(f"Error: VEP executable not found: {args.vep_path}")
        sys.exit(1)
    
    # Check if cache directory exists
    if not os.path.exists(args.cache_dir):
        print(f"Error: VEP cache directory not found: {args.cache_dir}")
        sys.exit(1)
    
    print("=== VEP Annotation and Filtering Pipeline ===")
    print(f"Input VCF: {args.input_vcf}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"VEP path: {args.vep_path}")
    print(f"Cache directory: {args.cache_dir}")
    print(f"Assembly: {args.assembly}")
    print(f"Threads: {args.threads}")
    print()
    
    try:
        # Step 1: Run VEP annotation
        print("Step 1: Running VEP annotation...")
        vep_output = run_vep_annotation(
            args.vep_path, 
            args.cache_dir, 
            args.input_vcf, 
            args.output_prefix,
            args.assembly,
            args.threads
        )
        
        # Step 2: Filter VEP output for target variants
        print("\nStep 2: Filtering for stop_gained and frameshift variants...")
        print("  - Keeping all variants with target consequences across all transcripts")
        
        filtered_output = f"{args.output_prefix}.stop_gained_frameshift.vcf"
        variant_count = filter_vep_output(vep_output, filtered_output)
        
        # Clean up full VEP output if not keeping it
        if not args.keep_full_vep:
            try:
                os.remove(vep_output)
                print(f"Removed full VEP output: {vep_output}")
            except OSError:
                print(f"Warning: Could not remove {vep_output}")
        
        print(f"\n=== Pipeline completed successfully! ===")
        print(f"Output files:")
        if args.keep_full_vep:
            print(f"  - Full VEP output: {vep_output}")
        print(f"  - Filtered VEP output (VCF): {filtered_output}")
        print(f"Found {variant_count} variants with stop_gained or frameshift consequences")
        
        if variant_count == 0:
            print("\nNote: No stop_gained or frameshift variants found in the input.")
        
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
