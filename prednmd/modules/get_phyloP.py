import pandas as pd
import argparse
import pyBigWig
import numpy as np

def clean_variant_data(df):
    """
    Clean and validate variant data types
    
    Args:
        df (pandas.DataFrame): Input dataframe
    
    Returns:
        pandas.DataFrame: Cleaned dataframe
    """
   
    clean_df = df.copy()
    
    # Convert chromosome to string
    clean_df['chrom'] = clean_df['CHR'].astype(str)
    # Remove 'chr' prefix if present
    clean_df['chrom'] = clean_df['chrom'].str.replace('chr', '')
    
    # Convert position to integer
    clean_df['pos'] = pd.to_numeric(clean_df['POS'], errors='coerce').astype('Int64')
    
    # Convert ref and alt to uppercase strings
    clean_df['ref'] = clean_df['REF_ALLELE'].astype(str).str.upper()
    clean_df['alt'] = clean_df['ALT_ALLELE'].astype(str).str.upper()
    
    return clean_df

def get_phylop_scores(variants_df, bigwig_path):
    """
    Get phyloP scores for variants using bigWig file
    
    Args:
        variants_df: pandas DataFrame with 'chrom', 'pos', 'ref', 'alt' columns
                    and optionally 'PTC_POS' for frameshift variants
        bigwig_path: Path to phyloP bigWig file
    """
    bw = pyBigWig.open(bigwig_path)
    scores = []
    snp_count = 0
    frameshift_count = 0
    
    for _, row in variants_df.iterrows():
        try:
            # Add 'chr' prefix if not present
            chrom = f"chr{row['chrom']}" if not str(row['chrom']).startswith('chr') else str(row['chrom'])
            
            # Determine if this is an SNP (both ref and alt are single nucleotides)
            is_snp = (len(str(row['ref'])) == 1 and len(str(row['alt'])) == 1)
            
            if is_snp:
                # Use variant position for SNPs
                score = bw.values(chrom, 
                                int(row['pos'])-1,  # 0-based coordinate
                                int(row['pos']))[0]
                snp_count += 1
            else:
                # Use PTC position for frameshift/indel variants
                if 'PTC_POS' in row and pd.notna(row['PTC_POS']):
                    ptc_pos = int(row['PTC_POS'])
                    # Get scores for the 3-nucleotide stop codon region
                    scores_3nt = bw.values(chrom, 
                                         ptc_pos - 1,  # 0-based coordinate
                                         ptc_pos + 2)  # +2 to include 3 nucleotides
                    # Calculate average, excluding NaN values
                    valid_scores = [s for s in scores_3nt if s is not None and not np.isnan(s)]
                    if valid_scores:
                        score = np.mean(valid_scores)
                    else:
                        score = np.nan
                    frameshift_count += 1
                else:
                    # No PTC_POS available for non-SNP variant
                    score = np.nan
                    frameshift_count += 1
            
            scores.append(score)
        except:
            scores.append(np.nan)
    
    bw.close()
    print(f"phyloP extraction: {snp_count} SNPs (variant position), {frameshift_count} frameshift/indels (PTC position)")
    return scores

def add_phylop_to_file(input_file, output_file, phylop_file):
    """
    Read input file and add phyloP information
    Automatically detects SNPs vs frameshift variants for phyloP extraction
    """
    df = pd.read_csv(input_file, sep='\t')
    print(f"Read {len(df)} variants from input file")
    
    clean_df = clean_variant_data(df)
    
    # Create query dataframe, always include PTC_POS if available
    query_df = pd.DataFrame({
        'chrom': clean_df['chrom'],
        'pos': clean_df['pos'],
        'ref': clean_df['ref'],
        'alt': clean_df['alt']
    })
    
    # Add PTC_POS to query_df if it exists in the original data
    if 'PTC_POS' in df.columns:
        query_df['PTC_POS'] = pd.to_numeric(df['PTC_POS'], errors='coerce').astype('Int64')
    
    # Get phyloP scores
    print("Getting phyloP scores (auto-detecting SNPs vs frameshift variants)...")
    df['phyloP'] = get_phylop_scores(query_df, phylop_file)
    df['has_phylop'] = df['phyloP'].notna().astype(int)
    
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Add phyloP information to variant file')
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Path to input txt file')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Path to output txt file')
    parser.add_argument('-p', '--phylop',
                        required=True,
                        help='Path to phyloP bigWig file')
    
    args = parser.parse_args()
    
    try:
        add_phylop_to_file(args.input, args.output, args.phylop)
            
    except FileNotFoundError:
        print(f"Error: Could not find input file {args.input}")
    except pd.errors.EmptyDataError:
        print(f"Error: Input file {args.input} is empty")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    main()