#!/usr/bin/env python3
"""
Command-line interface for DeepNMD pipeline
"""

import argparse
import sys
from pathlib import Path
from .modules.filter_vcf_by_gene import is_vep_annotated, filter_vcf_by_csq

from deepnmd import NMDPipeline, Config


def create_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="DeepNMD: Predict Nonsense-Mediated Decay (NMD) in genetic variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete pipeline from raw VCF (VEP annotation will be auto-detected and run if needed)
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml
  
  # Filter for specific gene 
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --gene BRCA1
  
  # Output separate feature table with all features and SHAP values
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --output-features
  
  # Add all features and SHAP values to VCF INFO field
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --full-vcf-annotation
  
  # Start from existing feature table
  deepnmd run -i features.txt -o output_dir -s sample1 -c config.yaml --from-features
  
  # Generate predictions only (no VCF output)
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --no-vcf-output
  
  # Skip PTC check for SNVs (but not frameshifts) in step 3
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --skip-ptc-check
  
  # Specify custom AF column for step 3
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --af-column gnomADg_AF
  

  
  # Generate template config file
  deepnmd init-config -o config.yaml

For detailed options, run: deepnmd run -h
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Run complete pipeline
    run_parser = subparsers.add_parser('run', help='Run NMD prediction pipeline')
    run_parser.add_argument('-i', '--input', required=True, 
                           help='Input file (VCF, VEP-annotated VCF, or feature table)')
    run_parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    run_parser.add_argument('-s', '--sample-name', required=True, help='Sample name')
    run_parser.add_argument('-c', '--config', help='Configuration YAML file')
    
    # Workflow options (human-readable!)
    workflow_group = run_parser.add_argument_group('Workflow options')
    workflow_group.add_argument('--skip-filtering', action='store_true',
                              help='Skip protein-coding region filtering')
    workflow_group.add_argument('--gene', type=str, metavar='GENE',
                              help='Filter variants to specific gene (by SYMBOL or Ensembl ID)')
    workflow_group.add_argument('--from-features', action='store_true',
                              help='Start from existing feature table (tab-delimited)')
    workflow_group.add_argument('--no-vcf-output', action='store_true',
                              help='Generate predictions table only (skip adding NMD results to INFO field of the VCF)')
    workflow_group.add_argument('--only-vcf-annotation', action='store_true',
                              help='Just add NMD predictions to VCF (requires --predictions-file)')
    workflow_group.add_argument('--predictions-file', metavar='FILE',
                              help='Pre-computed predictions file (for --only-vcf-annotation)')
    
    # Output options
    output_group = run_parser.add_argument_group('Output options')
    output_group.add_argument('--output-features', action='store_true',
                            help='Output a separate feature table with all features and SHAP values (prediction table will only contain essential columns)')
    output_group.add_argument('--full-vcf-annotation', action='store_true',
                            help='Add all features and SHAP values (n_terminal_escape_contrib, c_terminal_escape_contrib, general_escape_contrib) to VCF INFO field in addition to standard NMD annotations')
    
    # Feature extraction options
    feature_group = run_parser.add_argument_group('Feature extraction options')
    feature_group.add_argument('--skip-ptc-check', action='store_true',
                             help='Skip PTC check for SNVs (NOT frameshifts) in step 3. Assumes all SNVs annotated as stop_gained create PTCs.')
    feature_group.add_argument('--af-column', type=str, metavar='COLUMN',
                             help='Specify which AF column to use in step 3 (default: auto-detect gnomAD_AF or gnomADg_AF)')
    
    # File retention options (human-readable!)
    files_group = run_parser.add_argument_group('File retention options')
    files_group.add_argument('--keep-all', action='store_true',
                           help='Keep all intermediate files')
    files_group.add_argument('--keep-filtered-vcf', action='store_true',
                           help='Keep protein-coding filtered VCF')
    
    # Reference files
    ref_group = run_parser.add_argument_group('Reference files (override config)')
    ref_group.add_argument('--gtf', metavar='GTF_FILE', help='GTF annotation file')
    ref_group.add_argument('--genome', metavar='GENOME_FILE', help='Genome FASTA file')
    ref_group.add_argument('--cds', metavar='CDS_FILE', help='CDS FASTA file')
    
    # Annotation files
    annot_group = run_parser.add_argument_group('Annotation files (override config)')
    annot_group.add_argument('--gnomad', metavar='GNOMAD_FILE', help='gnomAD constraint metrics file')
    annot_group.add_argument('--phylop', metavar='PHYLOP_FILE', help='phyloP bigWig file')
    annot_group.add_argument('--m6a', metavar='M6A_FILE', help='m6A annotation file')
    annot_group.add_argument('--expression', metavar='EXPRESSION_FILE', help='Gene expression file')
    
    # VEP options
    vep_group = run_parser.add_argument_group('VEP options (override config)')
    vep_group.add_argument('--vep-path', help='Path to VEP executable')
    vep_group.add_argument('--vep-cache', help='VEP cache directory')
    vep_group.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], 
                          help='Genome assembly')
    
    # Model
    model_group = run_parser.add_argument_group('Model options (override config)')
    model_group.add_argument('--model-dir', help='Directory containing Random Forest models')
    
    # Runtime options
    runtime_group = run_parser.add_argument_group('Runtime options')
    runtime_group.add_argument('--threads', type=int, help='Number of threads')
    runtime_group.add_argument('--no-canonical', action='store_true',
                             help='Include non-canonical transcripts')
    
    # Initialize config
    init_parser = subparsers.add_parser('init-config', 
                                       help='Generate template configuration file')
    init_parser.add_argument('-o', '--output', default='deepnmd_config.yaml',
                           help='Output configuration file')
    
    # Version
    parser.add_argument('--version', action='version', version='DeepNMD 1.0.0')
    
    return parser


def cmd_run(args):
    """Run pipeline with human-readable workflow options"""
    # Capture the original user command
    original_command = ' '.join(sys.argv)
    
    # Create or load config
    if args.config:
        config = Config(args.config)
    else:
        config = Config()
    
    # Override config with command-line arguments
    if args.gtf:
        config.set('reference', 'gtf_file', value=args.gtf)
    if args.genome:
        config.set('reference', 'genome_fasta', value=args.genome)
    if args.cds:
        config.set('reference', 'cds_fasta', value=args.cds)
    if args.gnomad:
        config.set('annotation', 'gnomad_file', value=args.gnomad)
    if args.phylop:
        config.set('annotation', 'phylop_bigwig', value=args.phylop)
    if args.m6a:
        config.set('annotation', 'm6a_file', value=args.m6a)
    if args.expression:
        config.set('annotation', 'expression_file', value=args.expression)
    if args.vep_path:
        config.set('vep', 'vep_path', value=args.vep_path)
    if args.vep_cache:
        config.set('vep', 'cache_dir', value=args.vep_cache)
    if args.assembly:
        config.set('vep', 'assembly', value=args.assembly)
    if args.model_dir:
        config.set('model', 'model_dir', value=args.model_dir)
    if args.threads:
        config.set('runtime', 'threads', value=args.threads)
    if args.no_canonical:
        config.set('runtime', 'canonical_only', value=False)
    
    # Determine workflow based on human-readable options
    start_step = 1
    end_step = 8
    gene_filter = None
    
    # Handle workflow options - CHECK IN PRIORITY ORDER
    if args.from_features:
        start_step = 7  # Start from RF model (feature table is input)
        print(f"INFO: Starting from feature table (Step 7)")
    elif args.gene:
        print(f"INFO: Gene filtering enabled for {args.gene}")
        # Auto-detect VCF type using the updated detection logic
        is_vep = is_vep_annotated(args.input)
        
        if is_vep:
            # VEP-annotated: Use CSQ field
            print("INFO: Detected VEP-annotated file")
            print("INFO: Using CSQ/annotation field for filtering")
            
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            filtered_vcf = output_dir / f"filtered_{args.gene}.vcf"
            
            try:
                filter_vcf_by_csq(args.input, args.gene, str(filtered_vcf))
                args.input = str(filtered_vcf)
                start_step = 3  # Skip to feature extraction
                print(f"INFO: Starting from Step 3")
            except Exception as e:
                print(f"ERROR: {e}")
                return 1
        else:
            # Raw VCF: Use GTF in Step 1
            print("INFO: Detected raw VCF (not VEP-annotated)")
            print("INFO: Using GTF coordinates in Step 1")
            start_step = 1
            gene_filter = args.gene
    elif args.skip_filtering:
        start_step = 2  # Start from VEP
        print(f"INFO: Skipping protein-coding filter (Step 2)")
    else:
        # Auto-detect if VEP annotation is present
        # If VEP-annotated, skip VEP step (step 2)
        # Pipeline will handle this in the run() method
        pass
    
    if args.only_vcf_annotation:
        if not args.predictions_file:
            print("Error: --only-vcf-annotation requires --predictions-file", file=sys.stderr)
            return 1
        start_step = 8
        end_step = 8
    elif args.no_vcf_output:
        end_step = 7  # Stop after predictions
    
    # Determine which files to keep
    keep_files = [2]  # Always keep VEP annotated VCF
    if args.keep_all:
        keep_intermediate = True
    else:
        keep_intermediate = False
        if args.keep_filtered_vcf:
            keep_files.append(1)
    
    # Validate conflicting options
    if args.from_features and args.skip_filtering:
        print("Error: --from-features cannot be combined with --skip-filtering", 
              file=sys.stderr)
        return 1
    
    if args.gene and args.from_features:
        print("Error: --gene requires VCF input, cannot be used with --from-features",
              file=sys.stderr)
        return 1
    
    # Create and run pipeline
    try:
        pipeline = NMDPipeline(config, args.output_dir, args.sample_name, original_command=original_command)
        
        # Run with determined parameters
        outputs = pipeline.run(
            args.input,
            start_step=start_step,
            end_step=end_step,
            keep_intermediate=keep_intermediate,
            keep_files=keep_files,
            predictions_file=args.predictions_file if args.only_vcf_annotation else None,
            gene_filter=gene_filter,
            output_features=args.output_features if hasattr(args, 'output_features') else False,
            full_vcf_annotation=args.full_vcf_annotation if hasattr(args, 'full_vcf_annotation') else False,
            skip_ptc_check=args.skip_ptc_check if hasattr(args, 'skip_ptc_check') else False,
            af_column=args.af_column if hasattr(args, 'af_column') else None
        )
        
        print("\n" + "=" * 70)
        print("Pipeline completed successfully!")
        print("=" * 70)
        
        if 'final_vcf' in outputs:
            print(f"\nFinal annotated VCF: {outputs['final_vcf']}")
        if 'predictions_txt' in outputs:
            print(f"Predictions table: {outputs['predictions_txt']}")
        if 'features_txt_output' in outputs:
            print(f"Features table: {outputs['features_txt_output']}")
        
        # Show kept intermediate files if any
        if keep_intermediate or keep_files:
            print("\nIntermediate files retained:")
            for key, path in outputs.items():
                if key not in ['final_vcf', 'predictions_txt', 'features_txt_output']:
                    if Path(path).exists():
                        print(f"  - {path}")
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


def cmd_init_config(args):
    """Generate template configuration file"""
    template = """# DeepNMD Configuration File
# Edit this file with your specific paths and settings

# Ensembl reference files (REQUIRED)
reference:
  gtf_file: /path/to/reference.gtf
  genome_fasta: /path/to/genome.fa
  cds_fasta: /path/to/cds.fa

# Annotation databases (REQUIRED)
annotation:
  gnomad_file: /path/to/gnomad.constraint_metrics.tsv
  phylop_bigwig: /path/to/phyloP.bw
  m6a_file: /path/to/m6A_annotations.txt
  expression_file: /path/to/gene_expression.csv

# VEP configuration (REQUIRED if using VEP)
vep:
  vep_path: /path/to/vep #path to vep executable
  cache_dir: /path/to/.vep #path to vep cache
  assembly: GRCh37  # or GRCh38

# Machine learning model (REQUIRED)
model:
  model_dir: /path/to/RF_models/

# Runtime settings
runtime:
  threads: 32
  canonical_only: true
"""
    
    with open(args.output, 'w') as f:
        f.write(template)
    
    print(f"Template configuration file created: {args.output}")
    print("Please edit this file with your specific paths and settings.")
    return 0


def main():
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return 1
    
    # Route to appropriate command
    if args.command == 'run':
        return cmd_run(args)
    elif args.command == 'init-config':
        return cmd_init_config(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
