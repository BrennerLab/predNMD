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
  # Complete pipeline from raw VCF
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml
  
  # Filter for specific gene (skips Step 1, faster!)
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --gene BRCA1
  
  # If you already have VEP annotations
  deepnmd run -i vep_annotated.vcf -o output_dir -s sample1 -c config.yaml --skip-vep
  
  # Start from existing feature table
  deepnmd run -i features.txt -o output_dir -s sample1 -c config.yaml --from-features
  
  # Generate predictions only (no VCF output)
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --no-vcf-output
  
  # Keep specific intermediate files
  deepnmd run -i input.vcf.gz -o output_dir -s sample1 -c config.yaml --keep-features
  
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
                              help='Skip protein-coding region filtering (input already filtered)')
    workflow_group.add_argument('--gene', type=str, metavar='GENE',
                              help='Filter variants to specific gene (by SYMBOL or Ensembl ID) - skips Step 1')
    workflow_group.add_argument('--skip-vep', action='store_true',
                              help='Skip VCF filtering and VEP annotation (input already VEP-annotated)')
    workflow_group.add_argument('--from-features', action='store_true',
                              help='Start from existing feature table (input is features.txt)')
    workflow_group.add_argument('--no-vcf-output', action='store_true',
                              help='Generate predictions table only (skip VCF annotation)')
    workflow_group.add_argument('--only-vcf-annotation', action='store_true',
                              help='Only add predictions to VCF (requires --predictions-file)')
    workflow_group.add_argument('--predictions-file', metavar='FILE',
                              help='Pre-computed predictions file (for --only-vcf-annotation)')
    
    # File retention options (human-readable!)
    files_group = run_parser.add_argument_group('File retention options')
    files_group.add_argument('--keep-all', action='store_true',
                           help='Keep all intermediate files')
    files_group.add_argument('--keep-filtered-vcf', action='store_true',
                           help='Keep protein-coding filtered VCF')
    #files_group.add_argument('--keep-vep-vcf', action='store_true',
    #                       help='Keep VEP-annotated VCF')
    files_group.add_argument('--keep-features', action='store_true',
                           help='Keep feature extraction outputs')
    
    # Reference files
    ref_group = run_parser.add_argument_group('Reference files (override config)')
    ref_group.add_argument('--gtf', help='GTF annotation file')
    ref_group.add_argument('--genome', help='Genome FASTA file')
    ref_group.add_argument('--cds', help='CDS FASTA file')
    
    # Annotation files
    annot_group = run_parser.add_argument_group('Annotation files (override config)')
    annot_group.add_argument('--gnomad', help='gnomAD constraint metrics file')
    annot_group.add_argument('--phylop', help='phyloP bigWig file')
    annot_group.add_argument('--m6a', help='m6A annotation file')
    annot_group.add_argument('--expression', help='Gene expression file')
    
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
    parser.add_argument('--version', action='version', version='DeepNMD 1.1.0')
    
    return parser


def cmd_run(args):
    """Run pipeline with human-readable workflow options"""
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
        # Auto-detect VCF type
        is_vep = is_vep_annotated(args.input)
        
        if args.skip_vep or is_vep:
            # VEP-annotated: Use CSQ field
            print("INFO: Detected VEP-annotated VCF")
            print("INFO: Using CSQ field for filtering")
            
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
            print("INFO: Detected raw VCF")
            print("INFO: Using GTF coordinates in Step 1")
            start_step = 1
            gene_filter = args.gene
    elif args.skip_vep:
        start_step = 3  # Start from feature extraction
        print(f"INFO: Skipping VEP annotation (Step 3)")
    elif args.skip_filtering:
        start_step = 2  # Start from VEP
        print(f"INFO: Skipping protein-coding filter (Step 2)")
    
    
    if args.only_vcf_annotation:
        if not args.predictions_file:
            print("Error: --only-vcf-annotation requires --predictions-file", file=sys.stderr)
            return 1
        start_step = 8
        end_step = 8
    elif args.no_vcf_output:
        end_step = 7  # Stop after predictions
    
    # Determine which files to keep
    keep_files = [2]
    if args.keep_all:
        keep_intermediate = True
    else:
        keep_intermediate = False
        if args.keep_filtered_vcf:
            keep_files.append(1)
        #if args.keep_vep_vcf:
        #    keep_files.append(2)
        if args.keep_features:
            keep_files.extend([3, 6])  # Keep both feature outputs
    
    # Validate conflicting options
    if args.from_features and (args.skip_vep or args.skip_filtering):
        print("Error: --from-features cannot be combined with --skip-vep or --skip-filtering", 
              file=sys.stderr)
        return 1
    
    if args.gene and args.from_features:
        print("Error: --gene requires VCF input, cannot be used with --from-features",
              file=sys.stderr)
        return 1
    
    # Create and run pipeline
    try:
        pipeline = NMDPipeline(config, args.output_dir, args.sample_name)
        
        # Run with determined parameters
        outputs = pipeline.run(
            args.input,
            start_step=start_step,
            end_step=end_step,
            keep_intermediate=keep_intermediate,
            keep_files=keep_files,
            predictions_file=args.predictions_file if args.only_vcf_annotation else None,
            #gene_filter=args.gene if hasattr(args, 'gene') and args.gene else None
            gene_filter=gene_filter
        )
        
        print("\n" + "=" * 70)
        print("Pipeline completed successfully!")
        print("=" * 70)
        
        if 'final_vcf' in outputs:
            print(f"\nFinal annotated VCF: {outputs['final_vcf']}")
        if 'predictions_txt' in outputs:
            print(f"Predictions table: {outputs['predictions_txt']}")
        
        # Show kept intermediate files if any
        if keep_intermediate or keep_files:
            print("\nIntermediate files retained:")
            for key, path in outputs.items():
                if key not in ['final_vcf', 'predictions_txt']:
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

# Reference genome files (REQUIRED)
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
