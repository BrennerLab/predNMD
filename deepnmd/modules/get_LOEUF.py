import pandas as pd
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Add LOEUF scores to variant data.')
    parser.add_argument('--gnomad', default='/n/electric/data/yaqisu/gnomad.v4.1.constraint_metrics.tsv', 
                        help='Path to gnomAD constraint metrics file (default: /n/electric/data/yaqisu/gnomad.v4.1.constraint_metrics.tsv) (GRCh38)')
    parser.add_argument('--input', required=True, help='Path to input variant file')
    parser.add_argument('--output', required=True, help='Path to output file')
    args = parser.parse_args()

    # Step 1: Load the gnomAD constraint metrics file
    print(f"Loading gnomAD constraint metrics file from {args.gnomad}...")
    gnomad_df = pd.read_csv(args.gnomad, sep='\t')
    
    # Step 2: Create mapping for transcript ID to LOEUF
    print("Creating transcript to LOEUF mapping...")
    transcript_to_loeuf = {}
    
    for _, row in gnomad_df.iterrows():
        transcript_id = row['transcript']
        
        # Some cleaning to ensure we capture the transcript ID correctly
        if isinstance(transcript_id, str):
            transcript_id = transcript_id.strip()
            
            # Extract LOEUF (loss-of-function observed/expected upper bound fraction)
            if 'lof.oe_ci.upper' in row and not pd.isna(row['lof.oe_ci.upper']):
                transcript_to_loeuf[transcript_id] = row['lof.oe_ci.upper']
    
    print(f"Found {len(transcript_to_loeuf)} transcript IDs with LOEUF values.")
    
    # Step 3: Load the input dataset
    print(f"Loading input dataset from {args.input}...")
    input_df = pd.read_csv(args.input, sep='\t')
    
    # Step 4: Add LOEUF column to input dataset
    print("Adding LOEUF column to dataset...")
    
    # Identify the column containing transcript IDs (assuming 'Feature' based on original script)
    transcript_id_column = 'Feature'
    if transcript_id_column not in input_df.columns:
        # Try to find a column that might contain transcript IDs
        possible_columns = [col for col in input_df.columns if 'transcript' in col.lower() or 'feature' in col.lower()]
        if possible_columns:
            transcript_id_column = possible_columns[0]
            print(f"Using '{transcript_id_column}' as the transcript ID column")
        else:
            raise ValueError("Could not identify a column containing transcript IDs. Please specify the column name.")
    
    # Add LOEUF values
    input_df['LOEUF'] = input_df[transcript_id_column].apply(
        lambda feature_id: transcript_to_loeuf.get(feature_id, None)
    )
    
    # Step 5: Save the updated dataset
    print(f"Saving updated dataset to {args.output}...")
    input_df.to_csv(args.output, sep='\t', index=False)
    
    # Report statistics
    loeuf_count = input_df['LOEUF'].notna().sum()
    total_rows = len(input_df)
    
    print(f"Successfully added LOEUF values for {loeuf_count}/{total_rows} entries ({loeuf_count/total_rows:.1%})")


if __name__ == "__main__":
    main()
