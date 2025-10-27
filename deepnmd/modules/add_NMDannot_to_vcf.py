#!/usr/bin/env python3
"""
Add NMD annotations (nmd_trigger_probability, c_terminal_probability, n_terminal_probability) 
from a tab-delimited file into a VEP-annotated VCF file at the transcript level.

The script matches annotations by chromosome, position, reference allele, 
alternate allele, and transcript ID, then appends the probability fields to each 
transcript's CSQ entry.

Version: 5.2 - Updated to include three probability fields: NMD_PROB, C_TERMINAL_PROB, N_TERMINAL_PROB
"""

import sys
import argparse

def normalize_chrom(chrom):
    """Normalize chromosome name to include 'chr' prefix."""
    if not chrom.startswith('chr'):
        return 'chr' + chrom
    return chrom

def normalize_transcript_id(transcript_id):
    """Normalize transcript ID by removing version number."""
    if not transcript_id:
        return transcript_id
    # Remove version number (e.g., ENST00000600779.1 -> ENST00000600779)
    return transcript_id.split('.')[0]

def parse_annotation_file(annotation_file, debug=False):
    """
    Parse the tab-delimited annotation file and create a dictionary
    keyed by (chr, pos, ref, alt, transcript_id) tuple.
    """
    annotations = {}
    
    with open(annotation_file, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')
        
        if debug:
            print(f"\nAnnotation file has {len(header)} columns")
            print(f"Header: {header[:10]}...")  # Show first 10 columns
            print(f"Transcript IDs will be normalized (version numbers removed)")
        
        # Find column indices - try different possible names
        try:
            chr_idx = header.index('CHR')
            pos_idx = header.index('POS')
            ref_idx = header.index('REF_ALLELE')
            alt_idx = header.index('ALT_ALLELE')
            
            # Try to find transcript column - could be named differently
            transcript_idx = None
            for possible_name in ['transcript_id', 'Feature', 'TRANSCRIPT', 'Transcript']:
                if possible_name in header:
                    transcript_idx = header.index(possible_name)
                    if debug:
                        print(f"Found transcript column: '{possible_name}' at index {transcript_idx}")
                    break
            
            if transcript_idx is None:
                print("ERROR: Could not find transcript ID column.", file=sys.stderr)
                print(f"Available columns: {', '.join(header)}", file=sys.stderr)
                sys.exit(1)
            
            nmd_prob_idx = header.index('nmd_trigger_probability')
            c_term_prob_idx = header.index('c_terminal_probability')
            n_term_prob_idx = header.index('n_terminal_probability')
            
        except ValueError as e:
            print(f"Error: Required column not found in annotation file: {e}", file=sys.stderr)
            print(f"Available columns: {', '.join(header)}", file=sys.stderr)
            sys.exit(1)
        
        # Read data lines
        for line_num, line in enumerate(f, start=2):
            if line.strip():
                fields = line.strip().split('\t')
                
                if len(fields) <= max(chr_idx, pos_idx, ref_idx, alt_idx, transcript_idx):
                    if debug:
                        print(f"Warning: Line {line_num} has insufficient columns, skipping")
                    continue
                
                # Create key (normalize chromosome format and transcript ID)
                chrom = normalize_chrom(fields[chr_idx].strip())
                pos = fields[pos_idx].strip()
                ref = fields[ref_idx].strip()
                alt = fields[alt_idx].strip()
                transcript = normalize_transcript_id(fields[transcript_idx].strip())
                
                # Skip if essential fields are empty
                if not transcript or not pos:
                    if debug:
                        print(f"Warning: Line {line_num} has empty transcript or position, skipping")
                    continue
                
                key = (chrom, pos, ref, alt, transcript)
                
                # Get annotation values (handle empty values)
                nmd_prob = fields[nmd_prob_idx].strip() if nmd_prob_idx < len(fields) else ''
                c_term_prob = fields[c_term_prob_idx].strip() if c_term_prob_idx < len(fields) else ''
                n_term_prob = fields[n_term_prob_idx].strip() if n_term_prob_idx < len(fields) else ''
                
                annotations[key] = {
                    'nmd_trigger_probability': nmd_prob if nmd_prob else '.',
                    'c_terminal_probability': c_term_prob if c_term_prob else '.',
                    'n_terminal_probability': n_term_prob if n_term_prob else '.'
                }
    
    return annotations

def parse_csq_format(vcf_file):
    """
    Parse the CSQ format string from VCF header to get field positions.
    """
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                # Extract the format string
                format_start = line.find('Format: ') + 8
                format_end = line.find('"', format_start)
                format_string = line[format_start:format_end]
                csq_fields = format_string.split('|')
                
                # Find Feature index (transcript ID)
                try:
                    feature_idx = csq_fields.index('Feature')
                    return csq_fields, feature_idx
                except ValueError:
                    print("Warning: 'Feature' field not found in CSQ format", file=sys.stderr)
                    return csq_fields, None
            elif line.startswith('#CHROM'):
                break
    
    return None, None

def add_annotations_to_vcf(vcf_file, annotations, output_file, debug=False):
    """
    Read VCF file, add annotations to CSQ field for each transcript, and write output.
    Returns the number of unique annotations matched.
    """
    # Parse CSQ format to find transcript position
    csq_fields, feature_idx = parse_csq_format(vcf_file)
    
    if csq_fields is None:
        print("Error: Could not parse CSQ format from VCF header", file=sys.stderr)
        sys.exit(1)
    
    if debug:
        print(f"\nCSQ format has {len(csq_fields)} fields")
        print(f"Feature (transcript) is at index: {feature_idx}")
        print(f"Transcript IDs will be normalized (version numbers removed)")
        print(f"Using CSQ allele field ONLY for matching")
    
    matches_found = 0
    matched_annotation_keys = set()  # Track unique annotations matched
    variants_processed = 0
    
    # Statistics counters
    match_stats = {
        'direct': 0,
        'empty_alt': 0
    }
    
    with open(vcf_file, 'r') as infile, open(output_file, 'w') as outfile:
        csq_header_updated = False
        
        for line in infile:
            # Handle header lines
            if line.startswith('##INFO=<ID=CSQ'):
                # Update CSQ header to include new fields
                if not csq_header_updated:
                    format_start = line.find('Format: ') + 8
                    format_end = line.find('"', format_start)
                    old_format = line[format_start:format_end]
                    new_format = old_format + '|NMD_PROB|C_TERMINAL_PROB|N_TERMINAL_PROB'
                    line = line[:format_start] + new_format + line[format_end:]
                    csq_header_updated = True
                outfile.write(line)
                
            elif line.startswith('##') or line.startswith('#CHROM'):
                outfile.write(line)
                
            else:
                # Process data lines
                fields = line.strip().split('\t')
                
                if len(fields) < 8:
                    outfile.write(line)
                    continue
                
                # Normalize chromosome format for matching
                chrom_normalized = normalize_chrom(fields[0])
                pos = fields[1]
                ref = fields[3]
                info = fields[7]
                
                variants_processed += 1
                
                # Parse INFO field to find CSQ
                info_parts = info.split(';')
                new_info_parts = []
                
                for part in info_parts:
                    if part.startswith('CSQ='):
                        csq_value = part[4:]  # Remove 'CSQ='
                        csq_transcripts = csq_value.split(',')
                        
                        # Process each transcript
                        updated_transcripts = []
                        for transcript_csq in csq_transcripts:
                            csq_values = transcript_csq.split('|')
                            
                            # Get the allele from CSQ (field 0) - THE definitive source
                            csq_allele = csq_values[0] if len(csq_values) > 0 else ''
                            
                            # VEP uses '-' for deletions, normalize to empty string
                            if csq_allele == '-':
                                csq_allele = ''
                            
                            # Get transcript ID and normalize it
                            transcript_id = ''
                            if feature_idx is not None and feature_idx < len(csq_values):
                                transcript_id = normalize_transcript_id(csq_values[feature_idx])
                            
                            if debug and variants_processed <= 2:
                                print(f"\nVariant {variants_processed}: CSQ allele='{csq_allele}', transcript={transcript_id}")
                            
                            # Try to find a match
                            matched = False
                            nmd_prob = '.'
                            c_term_prob = '.'
                            n_term_prob = '.'
                            actual_key = None
                            
                            # Try 1: Direct match with CSQ allele
                            key = (chrom_normalized, pos, ref, csq_allele, transcript_id)
                            
                            if key in annotations:
                                nmd_prob = annotations[key]['nmd_trigger_probability']
                                c_term_prob = annotations[key]['c_terminal_probability']
                                n_term_prob = annotations[key]['n_terminal_probability']
                                actual_key = key
                                matched = True
                                match_stats['direct'] += 1
                                
                                if debug and matches_found < 5:
                                    print(f"  ✓ MATCH FOUND (direct)! NMD_PROB={nmd_prob}, C_TERMINAL_PROB={c_term_prob}, N_TERMINAL_PROB={n_term_prob}")
                            
                            # Try 2: Deletion - annotation has empty ALT
                            # VCF/CSQ: ACACT -> A (or -), Annotation: ACACT -> ""
                            elif len(ref) > 1 and len(csq_allele) == 1 and ref.startswith(csq_allele):
                                key_del = (chrom_normalized, pos, ref, '', transcript_id)
                                
                                if key_del in annotations:
                                    nmd_prob = annotations[key_del]['nmd_trigger_probability']
                                    c_term_prob = annotations[key_del]['c_terminal_probability']
                                    n_term_prob = annotations[key_del]['n_terminal_probability']
                                    actual_key = key_del
                                    matched = True
                                    match_stats['empty_alt'] += 1
                                    
                                    if debug and matches_found < 5:
                                        print(f"  ✓ MATCH FOUND (empty ALT)! NMD_PROB={nmd_prob}, C_TERMINAL_PROB={c_term_prob}, N_TERMINAL_PROB={n_term_prob}")
                            
                            # Track unique annotations
                            if matched and actual_key:
                                matched_annotation_keys.add(actual_key)
                                matches_found += 1
                            elif not matched and debug and variants_processed <= 2:
                                print(f"  ✗ No match found")
                            
                            # Append new fields to CSQ entry
                            csq_values.append(nmd_prob)
                            csq_values.append(c_term_prob)
                            csq_values.append(n_term_prob)
                            updated_transcripts.append('|'.join(csq_values))
                        
                        # Reconstruct CSQ field
                        new_info_parts.append('CSQ=' + ','.join(updated_transcripts))
                    else:
                        new_info_parts.append(part)
                
                fields[7] = ';'.join(new_info_parts)
                
                # Write the line
                outfile.write('\t'.join(fields) + '\n')
    
    if debug:
        print(f"\nProcessed {variants_processed} variants")
        print(f"Total match instances: {matches_found}")
        print(f"Unique annotations matched: {len(matched_annotation_keys)}")
    
    return len(matched_annotation_keys), match_stats  # Return unique count and stats

def main():
    parser = argparse.ArgumentParser(
        description='Add NMD annotations to VEP-annotated VCF file at the transcript level'
    )
    parser.add_argument(
        '-v', '--vcf',
        required=True,
        help='Input VEP-annotated VCF file'
    )
    parser.add_argument(
        '-a', '--annotations',
        required=True,
        help='Input tab-delimited annotation file with NMD data (must include transcript_id column)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output VCF file with added annotations'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Print debugging information'
    )
    
    args = parser.parse_args()
    
    print("Reading annotation file...")
    annotations = parse_annotation_file(args.annotations, args.debug)
    print(f"Loaded {len(annotations)} transcript-level annotations")
    
    if args.debug and len(annotations) > 0:
        print("\nFirst 5 annotation keys:")
        for i, (key, val) in enumerate(list(annotations.items())[:5]):
            print(f"  {key}")
            print(f"    NMD_PROB: {val['nmd_trigger_probability']}")
            print(f"    C_TERMINAL_PROB: {val['c_terminal_probability']}")
            print(f"    N_TERMINAL_PROB: {val['n_terminal_probability']}")
    
    print("\nProcessing VCF file and adding annotations to CSQ field...")
    matches_found, match_stats = add_annotations_to_vcf(args.vcf, annotations, args.output, args.debug)
    print(f"Found {matches_found} matching transcript annotations")
    
    # Print matching statistics
    print(f"\nMatching Strategy Statistics:")
    print(f"  Direct matches:        {match_stats['direct']} ({100*match_stats['direct']/matches_found:.1f}%)" if matches_found > 0 else "  Direct matches:        0")
    print(f"  Empty ALT matches:     {match_stats['empty_alt']} ({100*match_stats['empty_alt']/matches_found:.1f}%)" if matches_found > 0 else "  Empty ALT matches:     0")
    print(f"  Total:                 {sum(match_stats.values())}")
    
    print(f"\nDone! Output written to {args.output}")

if __name__ == '__main__':
    main()
