#!/usr/bin/env python3
"""
Gene filtering module for deepNMD

"""

import gzip
import sys
from pathlib import Path


def is_vep_annotated(vcf_file):
    """
    Check if VCF has VEP CSQ annotations.
    
    Args:
        vcf_file: Path to VCF file
        
    Returns:
        True if VCF has CSQ header, False otherwise
    """
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    try:
        with open_func(vcf_file, mode) as f:
            for line in f:
                if line.startswith('##INFO=<ID=CSQ'):
                    return True
                if line.startswith('#CHROM'):
                    return False
    except Exception as e:
        print(f"Warning: Could not check VCF annotation status: {e}")
        return False
    
    return False


def parse_csq_format(info_line):
    """
    Parse VEP CSQ format string to get field positions.
    
    Args:
        info_line: CSQ INFO header line
        
    Returns:
        Dictionary mapping field names to positions
    """
    if 'Format:' not in info_line:
        raise ValueError("Cannot find Format specification in CSQ header")
    
    format_part = info_line.split('Format:')[1].split('"')[0].strip()
    fields = [f.strip() for f in format_part.split('|')]
    
    return {field: idx for idx, field in enumerate(fields)}


def filter_vcf_by_csq(input_vcf, gene_name, output_vcf):
    """
    Filter VEP-annotated VCF using CSQ field.
    Fast filtering that preserves all VEP annotations.
    
    Args:
        input_vcf: Path to VEP-annotated VCF
        gene_name: Gene symbol or Ensembl ID (e.g., 'BRCA1' or 'ENSG00000012048')
        output_vcf: Path to output filtered VCF
        
    Returns:
        Number of variants found
        
    Raises:
        ValueError: If gene not found or VCF not VEP-annotated
    """
    print(f"Filtering VEP-annotated VCF for gene: {gene_name}")
    print(f"Input: {input_vcf}")
    print(f"Output: {output_vcf}")
    
    open_func = gzip.open if input_vcf.endswith('.gz') else open
    mode = 'rt' if input_vcf.endswith('.gz') else 'r'
    
    field_index = None
    variants_found = 0
    has_csq_header = False
    
    with open_func(input_vcf, mode) as infile, \
         open(output_vcf, 'w') as outfile:
        
        for line in infile:
            # Handle header lines
            if line.startswith('##'):
                outfile.write(line)
                
                # Parse CSQ format
                if line.startswith('##INFO=<ID=CSQ'):
                    has_csq_header = True
                    try:
                        field_index = parse_csq_format(line)
                        print(f"Found VEP CSQ fields: {', '.join(field_index.keys())}")
                    except Exception as e:
                        raise ValueError(f"Could not parse CSQ format: {e}")
                
                continue
            
            # Column header line
            if line.startswith('#CHROM'):
                if not has_csq_header:
                    raise ValueError(
                        "VCF does not appear to be VEP-annotated. "
                        "Missing ##INFO=<ID=CSQ header line. "
                        "Run VEP first or use raw VCF filtering."
                    )
                outfile.write(line)
                continue
            
            # Variant lines
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            info = fields[7]
            
            # Extract CSQ annotation
            csq_value = None
            for info_field in info.split(';'):
                if info_field.startswith('CSQ='):
                    csq_value = info_field[4:]  # Remove 'CSQ='
                    break
            
            if not csq_value or csq_value == '.':
                continue
            
            # Check all transcript annotations for gene match
            matched = False
            transcripts = csq_value.split(',')
            
            for transcript in transcripts:
                csq_fields = transcript.split('|')
                
                # Check SYMBOL field
                if 'SYMBOL' in field_index:
                    idx = field_index['SYMBOL']
                    if idx < len(csq_fields):
                        symbol = csq_fields[idx]
                        if symbol == gene_name:
                            matched = True
                            break
                
                # Check Gene field (Ensembl ID)
                if 'Gene' in field_index:
                    idx = field_index['Gene']
                    if idx < len(csq_fields):
                        gene_id = csq_fields[idx]
                        if gene_id == gene_name:
                            matched = True
                            break
            
            if matched:
                outfile.write(line)
                variants_found += 1
    
    print(f"Filtered {variants_found} variants for gene {gene_name}")
    
    if variants_found == 0:
        raise ValueError(
            f"No variants found for gene '{gene_name}'.\n"
            f"Please check:\n"
            f"  1. Gene name spelling (case-sensitive)\n"
            f"  2. Try Ensembl ID instead (e.g., ENSG00000012048 for BRCA1)\n"
            f"  3. VCF contains variants in this gene:\n"
            f"     grep 'SYMBOL={gene_name}' {input_vcf}"
        )
    
    return variants_found


def filter_gene(input_vcf, gene_name, output_vcf, gtf_file=None, force_mode=None):
    """
    Smart gene filtering that auto-detects VCF type.
    
    Args:
        input_vcf: Input VCF file
        gene_name: Gene symbol or Ensembl ID
        output_vcf: Output filtered VCF
        gtf_file: GTF file (optional, only for raw VCF mode)
        force_mode: Force 'csq' or 'gtf' mode (optional, auto-detects if None)
        
    Returns:
        Tuple of (filter_mode_used, variants_found)
    """
    # Determine filter mode
    if force_mode:
        mode = force_mode
    elif is_vep_annotated(input_vcf):
        mode = 'csq'
    else:
        mode = 'gtf'
    
    print(f"Gene filtering mode: {mode.upper()}")
    
    if mode == 'csq':
        # Use CSQ field filtering 
        variants = filter_vcf_by_csq(input_vcf, gene_name, output_vcf)
        return ('csq', variants)
    
    elif mode == 'gtf':
        # Use GTF coordinate filtering
        # This should be handled by Step 1 script with --gene option
        raise NotImplementedError(
            "GTF mode should be handled by Step 1 script (vcf_filterProteinCoding.py) "
            "with --gene option. This function is for VEP-annotated VCFs only."
        )
    
    else:
        raise ValueError(f"Unknown filter mode: {mode}")


def main():
    """Command-line interface for gene filtering"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Filter VCF for specific gene (VEP-annotated VCFs only)',
        epilog="""
Note: For raw (non-VEP-annotated) VCFs, use vcf_filterProteinCoding.py with --gene option instead.
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input VEP-annotated VCF file')
    parser.add_argument('-g', '--gene', required=True,
                       help='Gene name (e.g., BRCA1) or Ensembl ID (e.g., ENSG00000012048)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output filtered VCF file')
    
    args = parser.parse_args()
    
    try:
        # Check if VCF is VEP-annotated
        if not is_vep_annotated(args.input):
            print("Error: Input VCF does not appear to be VEP-annotated.")
            print("For raw VCFs, use vcf_filterProteinCoding.py with --gene option instead.")
            sys.exit(1)
        
        # Filter VCF
        filter_vcf_by_csq(args.input, args.gene, args.output)
        
        print("\nFiltering completed successfully!")
        print(f"Output saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
