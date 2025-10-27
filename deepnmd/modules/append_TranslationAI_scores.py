#!/usr/bin/env python3

import csv
import sys
import re
import subprocess
import os

def run_translationai(fasta_file):
    """
    Run TranslationAI on the input FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file
    
    Returns:
        Path to the generated TranslationAI output file, or None on error
    """
    print(f"\nRunning TranslationAI on {fasta_file}...")
    print("=" * 60)
    
    # Construct the command
    cmd = ["translationai", "-I", fasta_file, "-t", "0.0001,0.0001"]
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        # Print stdout if available
        if result.stdout:
            print("TranslationAI output:")
            print(result.stdout)
        
        # Construct expected output filename
        base_name = os.path.splitext(fasta_file)[0]
        output_file = f"{base_name}_predORFs_0.0001_0.0001.txt"
        
        # Verify the output file was created
        if os.path.exists(output_file):
            print(f"✓ TranslationAI completed successfully!")
            print(f"✓ Output file created: {output_file}")
            return output_file
        else:
            print(f"Error: Expected output file not found: {output_file}")
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"Error running TranslationAI: {e}")
        if e.stderr:
            print(f"Error output: {e.stderr}")
        return None
    except FileNotFoundError:
        print("Error: 'translationai' command not found. Please ensure TranslationAI is installed and in your PATH.")
        return None
    except Exception as e:
        print(f"Unexpected error running TranslationAI: {e}")
        return None

def parse_translationai_output(translationai_file):
    """
    Parse TranslationAI output and extract TIS scores for downstream AUG and TTS scores for PTC.
    
    Args:
        translationai_file: Path to TranslationAI output file
    
    Returns:
        Dictionary mapping sequence identifiers to scores
    """
    print(f"\nParsing TranslationAI output from {translationai_file}...")
    scores_dict = {}
    
    try:
        with open(translationai_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Split by tab: identifier and score entries
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                
                seq_id = parts[0].strip()
                score_parts = parts[1:]
                
                # Extract downstream AUG and PTC positions from header
                # Format: >chr1:889455-889455(-)(ENST00000327044)(876, 705,)
                match = re.search(r'\((\d+|0), (\d+),\)', seq_id)
                if not match:
                    print(f"Warning: Cannot parse positions from {seq_id}")
                    continue
                
                downstream_aug_pos = int(match.group(1))
                ptc_pos = int(match.group(2))
                
                # Parse all score entries to find matching positions
                downstream_aug_tis_score = None
                ptc_tts_score = None
                
                # Combine all score parts and split by whitespace
                all_scores = ' '.join(score_parts).split()
                
                for entry in all_scores:
                    if ',' in entry:
                        try:
                            tis_pos, tts_pos, tis_score, tts_score = entry.split(',')
                            tis_pos = int(tis_pos)
                            tts_pos = int(tts_pos)
                            tis_score = float(tis_score)
                            tts_score = float(tts_score)
                            
                            # Match downstream AUG position for TIS score
                            if downstream_aug_pos > 0 and tis_pos == downstream_aug_pos:
                                downstream_aug_tis_score = tis_score
                            
                            # Match PTC position for TTS score
                            if tts_pos == ptc_pos:
                                ptc_tts_score = tts_score
                                
                        except (ValueError, IndexError):
                            continue
                
                # Store scores
                scores_dict[seq_id] = {
                    'downstream_aug_pos': downstream_aug_pos,
                    'ptc_pos': ptc_pos,
                    'downstream_aug_tis_score': downstream_aug_tis_score,
                    'ptc_tts_score': ptc_tts_score
                }
    
    except FileNotFoundError:
        print(f"Error: TranslationAI file '{translationai_file}' not found")
        return None
    except Exception as e:
        print(f"Error parsing TranslationAI file: {e}")
        return None
    
    print(f"✓ Parsed {len(scores_dict)} sequences with scores")
    return scores_dict

def create_sequence_identifier_pattern(row):
    """
    Create sequence identifier pattern to match FASTA header.
    
    Args:
        row: Dictionary containing row data
    
    Returns:
        Base pattern string for matching
    """
    try:
        chr_name = f"chr{row['CHR']}"
        pos = int(float(row['POS']))
        
        # Try both 'Strand' and 'strand' for flexibility
        strand = row.get('Strand', row.get('strand', ''))
        
        # Try multiple possible column names for transcript/feature ID
        feature = row.get('transcript_id', row.get('Feature', ''))
        
        # Create base pattern that should be in the FASTA header
        pattern = f">{chr_name}:{pos}-{pos}({strand})({feature})"
        return pattern
        
    except (ValueError, TypeError, KeyError) as e:
        print(f"Warning: Could not create pattern for row: {e}")
        return None

def append_scores_to_original(original_file, translationai_file, output_file):
    """
    Append TranslationAI scores to original file.
    
    Args:
        original_file: Original TXT file
        translationai_file: TranslationAI output file
        output_file: Output file with appended scores
    """
    # Parse TranslationAI scores
    scores_dict = parse_translationai_output(translationai_file)
    if scores_dict is None:
        return False
    
    print(f"\nAppending scores to {original_file}...")
    
    try:
        with open(original_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            reader = csv.DictReader(infile, delimiter='\t')
            
            # Add new columns
            fieldnames = list(reader.fieldnames) + [
                'downstream_inframeAUG_translationAI',
                'PTC_translationAI'
            ]
            
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            matched_count = 0
            total_count = 0
            
            for row in reader:
                total_count += 1
                
                # Create identifier pattern
                pattern = create_sequence_identifier_pattern(row)
                if pattern is None:
                    row['downstream_inframeAUG_translationAI'] = 0.0001
                    row['PTC_translationAI'] = 0.0001
                    writer.writerow(row)
                    continue
                
                # Find matching sequence in scores
                matched_seq_id = None
                for seq_id in scores_dict:
                    if pattern in seq_id:
                        matched_seq_id = seq_id
                        break
                
                if matched_seq_id:
                    scores = scores_dict[matched_seq_id]
                    matched_count += 1
                    
                    # Add downstream AUG TIS score (use 0.0001 for missing scores)
                    if scores['downstream_aug_tis_score'] is not None:
                        row['downstream_inframeAUG_translationAI'] = scores['downstream_aug_tis_score']
                    else:
                        row['downstream_inframeAUG_translationAI'] = 0.0001
                    
                    # Add PTC TTS score (use 0.0001 for missing scores)
                    if scores['ptc_tts_score'] is not None:
                        row['PTC_translationAI'] = scores['ptc_tts_score']
                    else:
                        row['PTC_translationAI'] = 0.0001
                else:
                    # No match found - use 0.0001 for both scores
                    row['downstream_inframeAUG_translationAI'] = 0.0001
                    row['PTC_translationAI'] = 0.0001
                
                writer.writerow(row)
            
            print(f"\n{'=' * 60}")
            print(f"✓ Score appending completed!")
            print(f"  Total rows processed: {total_count}")
            print(f"  Rows with matched scores: {matched_count}")
            print(f"  Match rate: {matched_count/total_count*100:.1f}%")
            print(f"  Output written to: {output_file}")
            print(f"{'=' * 60}")
            return True
            
    except FileNotFoundError:
        print(f"Error: Original file '{original_file}' not found")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False

def main():
    if len(sys.argv) != 4:
        print("TranslationAI Score Appender with Integrated Runner")
        print("=" * 60)
        print("Usage: python script.py <original.txt> <input.fasta> <output.txt>")
        print("\nThis script will:")
        print("  1. Run TranslationAI on the input FASTA file")
        print("  2. Parse the TranslationAI output")
        print("  3. Append scores to the original TXT file:")
        print("     - downstream_inframeAUG_translationAI (TIS score)")
        print("     - PTC_translationAI (TTS score)")
        print("\nExample:")
        print("  python script.py variants.txt sequences.fasta output_with_scores.txt")
        print("\nNote: TranslationAI must be installed and available in your PATH")
        sys.exit(1)
    
    original_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Verify input files exist
    if not os.path.exists(original_file):
        print(f"Error: Original file '{original_file}' not found")
        sys.exit(1)
    
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file '{fasta_file}' not found")
        sys.exit(1)
    
    print("TranslationAI Score Appender")
    print("=" * 60)
    print(f"Original file: {original_file}")
    print(f"FASTA file: {fasta_file}")
    print(f"Output file: {output_file}")
    
    # Step 1: Run TranslationAI
    translationai_output = run_translationai(fasta_file)
    if translationai_output is None:
        print("\nError: Failed to run TranslationAI")
        sys.exit(1)
    
    # Step 2: Append scores
    success = append_scores_to_original(original_file, translationai_output, output_file)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
