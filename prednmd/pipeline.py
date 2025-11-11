"""
Main pipeline orchestrator for predNMD
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any
import logging

from .config import Config


class NMDPipeline:
    """
    Complete NMD prediction pipeline orchestrator
    
    This class manages the entire workflow from VCF input to NMD predictions.
    """
    
    def __init__(self, config: Config, output_dir: str, sample_name: str, original_command: str = None):
        """
        Initialize pipeline
        
        Args:
            config: Configuration object
            output_dir: Output directory for all results
            sample_name: Sample name for output files
        """
        self.config = config
        self.output_dir = Path(output_dir)
        self.sample_name = sample_name
        self.original_command = original_command
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.logger = self._setup_logging()
        
        # Get module directory
        self.module_dir = Path(__file__).parent / 'modules'
    
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        logger = logging.getLogger('predNMD')
        logger.setLevel(logging.INFO)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # File handler
        log_file = self.output_dir / f"{self.sample_name}_prednmd.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)
        
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)
        
        return logger
    
    def _detect_input_type(self, input_file: str) -> str:
        """
        Detect the type of input file
        
        Returns:
            'raw_vcf', 'vep_vcf', 'vep_tab_extra', 'vep_tab_standard', or 'feature_table'
        """
        # Check if it's a feature table (TSV/TXT)
        if input_file.endswith('.txt') or input_file.endswith('.tsv'):
            # Check if it's actually a VEP tab-delimited file
            try:
                with open(input_file, 'r') as f:
                    first_line = f.readline()
                    # VEP tab files have specific headers
                    if '#Uploaded_variation' in first_line or 'Uploaded_variation' in first_line:
                        if '\tExtra' in first_line or first_line.rstrip().endswith('Extra'):
                            return 'vep_tab_extra'
                        else:
                            return 'vep_tab_standard'
                # Otherwise it's a feature table
                return 'feature_table'
            except:
                return 'feature_table'
       
        else:
            # It's a VCF file (could be .vcf, .vcf.gz, etc.)
            # Check for VEP annotation in VCF
            try:
                if input_file.endswith('.gz'):
                    import gzip
                    with gzip.open(input_file, 'rt') as f:
                        for line in f:
                            if line.startswith('##INFO=<ID=CSQ'):
                                return 'vep_vcf'
                            if not line.startswith('#'):
                                break
                else:
                    with open(input_file, 'r') as f:
                        for line in f:
                            if line.startswith('##INFO=<ID=CSQ'):
                                return 'vep_vcf'
                            if not line.startswith('#'):
                                break
                return 'raw_vcf'
            except:
                return 'unknown'
  
    
    def _validate_config_for_steps(self, start_step: int, end_step: int):
        """
        Validate that required configuration is present for the steps being run
        
        Args:
            start_step: First step to run
            end_step: Last step to run
        """
        required_by_step = {
            1: [('reference', 'gtf_file')],
            2: [('vep', 'cache_dir')],
            3: [('reference', 'gtf_file'), ('reference', 'genome_fasta'), #('reference', 'cds_fasta'),
                ('annotation', 'm6a_file'), ('annotation', 'expression_file')],
            4: [('annotation', 'gnomad_file')],
            5: [('annotation', 'phylop_bigwig')],
            6: [],  # TranslationAI has no config requirements
            7: [('model', 'model_dir')],
            8: [],  # VCF annotation has no config requirements
        }
        
        missing = []
        for step in range(start_step, end_step + 1):
            for keys in required_by_step.get(step, []):
                if self.config.get(*keys) is None:
                    missing.append(f"Step {step} requires: {'.'.join(keys)}")
        
        if missing:
            raise ValueError(f"Missing required configuration:\n  " + "\n  ".join(missing))
        
        return True
    
    def _run_step(self, script_name: str, args: list, step_name: str) -> str:
        """
        Run a pipeline step
        
        Args:
            script_name: Name of the script to run
            args: Arguments to pass to the script
            step_name: Name of the step for logging
        
        Returns:
            Output file path from the step
        """
        self.logger.info(f"Starting Step: {step_name}")
        self.logger.info(f"Command: python {script_name} {' '.join(map(str, args))}")
        
        script_path = self.module_dir / script_name
        cmd = ['python', str(script_path)] + [str(arg) for arg in args]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"Completed Step: {step_name}")
            if result.stdout:
                self.logger.debug(f"STDOUT: {result.stdout}")
            return args[-1] if args else None  # Usually last arg is output
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed Step: {step_name}")
            self.logger.error(f"Return code: {e.returncode}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise RuntimeError(f"Pipeline failed at step: {step_name}")
    
    def run(self, input_file: str, 
            start_step: int = 1,
            end_step: int = 8,
            keep_intermediate: bool = None,
            keep_files: Optional[list] = None,
            predictions_file: Optional[str] = None,
            gene_filter: Optional[str] = None,
            output_features: bool = False,
            full_vcf_annotation: bool = False,
            skip_ptc_check: bool = False,
            af_column: Optional[str] = None) -> Dict[str, str]:
        """
        Run the pipeline with flexible entry and exit points
        
        Args:
            input_file: Path to input file (VCF, annotated VCF, or feature table)
            start_step: Pipeline step to start from (1-8)
            end_step: Pipeline step to end at (1-8)
            keep_intermediate: Whether to keep all intermediate files (overrides config)
            keep_files: Specific intermediate files to keep (list of step numbers)
            predictions_file: Pre-computed predictions file (for step 8 only mode)
            gene_filter: Optional gene name to filter variants (SYMBOL or Ensembl ID)
            output_features: If True, output separate feature table with all features and SHAP values
            full_vcf_annotation: If True, add all features and SHAP values to VCF INFO field
            skip_ptc_check: If True, skip PTC check for SNVs (not frameshifts) in step 3
            af_column: Specific AF column to use in step 3 (default: auto-detect)
        
        Returns:
            Dictionary with paths to output files
        """
        if keep_intermediate is None:
            keep_intermediate = self.config.get('runtime', 'keep_intermediate', default=False)
        
        if keep_files is None:
            keep_files = []
        
        self.logger.info("=" * 70)
        self.logger.info("predNMD Pipeline Started")
        self.logger.info("=" * 70)
        self.logger.info(f"Input file: {input_file}")
        self.logger.info(f"Sample: {self.sample_name}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Running steps: {start_step} to {end_step}")
        
        # Validate config based on steps being run
        try:
            self._validate_config_for_steps(start_step, end_step)
        except ValueError as e:
            self.logger.error(f"Configuration validation failed: {e}")
            raise
        
        outputs = {}
        threads = self.config.get('runtime', 'threads', default=1)
        current_file = input_file

        # Detect input file type
        input_type = self._detect_input_type(input_file)
        self.logger.info(f"Detected input file type: {input_type}")
        
        # Auto-adjust start_step based on input type if starting from step 1 or 2
        if start_step <= 2:
            if input_type in ['vep_vcf', 'vep_tab_extra', 'vep_tab_standard']:
                # Input is already VEP-annotated, skip VEP annotation step
                if start_step == 1:
                    self.logger.info("Input is VEP-annotated. Skipping protein-coding filter (Step 1) and VEP annotation (Step 2).")
                    self.logger.info("Starting from Step 3 (feature extraction).")
                    start_step = 3
                elif start_step == 2:
                    self.logger.info("Input is VEP-annotated. Skipping VEP annotation (Step 2).")
                    self.logger.info("Starting from Step 3 (feature extraction).")
                    start_step = 3
            else:
                # Raw VCF - need to go through VEP annotation
                self.logger.info("Input is raw VCF. Will run VEP annotation in Step 2.")
        
        # Step 1: Filter VCF for protein-coding regions
        if start_step <= 1 <= end_step:
            protein_coding_vcf = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz"
            if gene_filter:
                self._run_step(
                    'vcf_filterProteinCoding.py',
                    [
                        current_file,
                        protein_coding_vcf,
                        '--gtf-file', self.config.get('reference', 'gtf_file'),
                        '--threads', threads,
                        '--gene', gene_filter
                    ],
                    "1. Filter protein-coding regions with gene filter"
                )
            else:
                self._run_step(
                    'vcf_filterProteinCoding.py',
                    [
                        current_file,
                        protein_coding_vcf,
                        '--gtf-file', self.config.get('reference', 'gtf_file'),
                        '--threads', threads
                    ],
                    "1. Filter protein-coding regions"
                )
            outputs['protein_coding_vcf'] = str(protein_coding_vcf)
            current_file = protein_coding_vcf
        
        # Step 2: VEP annotation
        if start_step <= 2 <= end_step:
            vep_prefix = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz"
            vep_args = [
                current_file,
                vep_prefix,
                '--vep-path', self.config.get('vep', 'vep_path', default='vep'),
                '--cache-dir', self.config.get('vep', 'cache_dir'),
                '--assembly', self.config.get('vep', 'assembly', default='GRCh37'),
                '--threads', threads
            ]
            if keep_intermediate or 2 in keep_files:
                vep_args.append('--keep-full-vep')
            
            self._run_step(
                'vepAnnot_filter.py',
                vep_args,
                "2. VEP annotation and filtering"
            )
            
            # Track both the complete VEP annotated VCF and the filtered VCF
            vep_annotated_vcf = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz.vep.vcf"
            stop_frameshift_vcf = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz.stop_gained_frameshift.vcf"
            outputs['vep_annotated_vcf'] = str(vep_annotated_vcf)
            outputs['vep_filtered_vcf'] = str(stop_frameshift_vcf)
            current_file = stop_frameshift_vcf
        
        # Gene filtering (if requested) - applies after VEP annotation
        if gene_filter and start_step <= 2:
            from .modules.filter_vcf_by_gene import filter_vcf_by_csq
            
            gene_filtered_vcf = self.output_dir / f"{self.sample_name}.{gene_filter}.vcf"
            self.logger.info(f"Filtering variants for gene: {gene_filter}")
            
            n_variants = filter_vcf_by_csq(str(current_file), gene_filter, str(gene_filtered_vcf))
            
            if n_variants == 0:
                raise ValueError(f"No variants found for gene '{gene_filter}'")
            
            self.logger.info(f"Filtered to {n_variants} variants in gene '{gene_filter}'")
            outputs['gene_filtered_vcf'] = str(gene_filtered_vcf)
            current_file = gene_filtered_vcf
        
        # Step 3: Add features
        if start_step <= 3 <= end_step:
            feature_txt = self.output_dir / f"{self.sample_name}.feature_added.txt"
            translationai_fasta = self.output_dir / f"{self.sample_name}.translationAI.fa"
            if self.config.get('reference', 'cds_fasta') is None:
                feature_args = [
                    current_file,
                    self.config.get('reference', 'gtf_file'),
                    self.config.get('reference', 'genome_fasta'),
                    '-o', feature_txt,
                    '--translationAI-fasta', translationai_fasta
                ]
            else:
                feature_args = [
                    current_file,
                    self.config.get('reference', 'gtf_file'),
                    self.config.get('reference', 'genome_fasta'),
                    '--cds-fasta', self.config.get('reference', 'cds_fasta'),
                    '-o', feature_txt,
                    '--translationAI-fasta', translationai_fasta
                ]
            
            # Optional annotations
            if self.config.get('annotation', 'm6a_file'):
                feature_args.extend(['--m6a-file', self.config.get('annotation', 'm6a_file')])
            if self.config.get('annotation', 'expression_file'):
                feature_args.extend(['--expression-file', self.config.get('annotation', 'expression_file')])
            if self.config.get('runtime', 'canonical_only', default=True):
                feature_args.append('--canonical')
            
            # Add skip-ptc-check option if specified
            if skip_ptc_check:
                feature_args.append('--skip-ptc-check')
                self.logger.info("PTC check will be skipped for SNVs (not frameshifts)")
            
            # Add AF column option if specified
            if af_column:
                feature_args.extend(['--af-col', af_column])
                self.logger.info(f"Using AF column: {af_column}")
            
            self._run_step(
                'check_ptc_add_features.py',
                feature_args,
                "3. Add features"
            )
            outputs['features_txt'] = str(feature_txt)
            outputs['translationai_fasta'] = str(translationai_fasta)
            current_file = feature_txt
        
        # Step 4: Add LOEUF
        if start_step <= 4 <= end_step:
            self._run_step(
                'get_LOEUF.py',
                [
                    '--gnomad', self.config.get('annotation', 'gnomad_file'),
                    '--input', current_file,
                    '--output', current_file
                ],
                "4. Add LOEUF annotation"
            )
        
        # Step 5: Add phyloP
        if start_step <= 5 <= end_step:
            self._run_step(
                'get_phyloP.py',
                [
                    '-i', current_file,
                    '-o', current_file,
                    '-p', self.config.get('annotation', 'phylop_bigwig')
                ],
                "5. Add phyloP annotation"
            )
        
        # Step 6: Add TranslationAI scores
        if start_step <= 6 <= end_step:
            feature_translationai_txt = self.output_dir / f"{self.sample_name}.feature_added_translationAI.txt"
            # Need the FASTA file from step 3
            if 'translationai_fasta' not in outputs:
                translationai_fasta = self.output_dir / f"{self.sample_name}.translationAI.fa"
                if not os.path.exists(translationai_fasta):
                    self.logger.error("TranslationAI FASTA file not found. Step 3 must be run first or file must exist.")
                    raise FileNotFoundError(f"Required file not found: {translationai_fasta}")
            else:
                translationai_fasta = outputs['translationai_fasta']
            
            self._run_step(
                'append_TranslationAI_scores.py',
                [
                    current_file,
                    translationai_fasta,
                    feature_translationai_txt
                ],
                "6. Add TranslationAI scores"
            )
            outputs['features_translationai_txt'] = str(feature_translationai_txt)
            current_file = feature_translationai_txt
        
        # Step 7: Apply Random Forest model
        if start_step <= 7 <= end_step:
            feature_pred_txt = self.output_dir / f"{self.sample_name}.with_predictions.txt"
            
            rf_args = [
                self.config.get('model', 'model_dir'),
                current_file,
                feature_pred_txt
            ]
            
            # Add --features-output option if requested
            if output_features:
                features_output_txt = self.output_dir / f"{self.sample_name}.features.txt"
                rf_args.extend(['--features-output', features_output_txt])
                self.logger.info("Separate features table will be generated")
            
            # Add original command if available
            if self.original_command:
                rf_args.extend(['--command', self.original_command])
            
            self._run_step(
                'get_RF_SHAP_N-escape.py',
                rf_args,
                "7. Apply Random Forest model"
            )
            outputs['predictions_txt'] = str(feature_pred_txt)
            if output_features:
                outputs['features_txt_output'] = str(features_output_txt)
            current_file = feature_pred_txt
        
        # Step 8: Add NMD annotation to VCF
        if start_step <= 8 <= end_step:
            # Determine which predictions file to use
            if predictions_file:
                # Use provided predictions file
                predictions_to_use = predictions_file
                self.logger.info(f"Using provided predictions file: {predictions_file}")
            elif 'predictions_txt' in outputs:
                # Use predictions from step 7
                predictions_to_use = outputs['predictions_txt']
            else:
                # Look for predictions file from previous run
                predictions_to_use = current_file
            
            # Use complete VEP annotated VCF (not the filtered one)
            # Priority: 1) VEP annotated from step 2, 2) User-provided VEP-annotated input, 3) Original input
            vcf_to_annotate = None
            
            if 'vep_annotated_vcf' in outputs:
                # Use complete VEP annotated VCF from step 2
                vcf_to_annotate = outputs['vep_annotated_vcf']
                self.logger.info(f"Using complete VEP annotated VCF from step 2: {vcf_to_annotate}")
            elif input_type in ['vep_vcf']:
                # User provided a VEP-annotated VCF directly
                vcf_to_annotate = input_file
                self.logger.info(f"Using user-provided VEP annotated VCF: {vcf_to_annotate}")
            else:
                self.logger.warning("Step 8 requires a VEP-annotated VCF file. Skipping VCF annotation.")
            
            if vcf_to_annotate:
                final_vcf = self.output_dir / f"{self.sample_name}.NMDannot.vcf"
                
                vcf_args = [
                    '-v', vcf_to_annotate,
                    '-a', predictions_to_use,
                    '-o', final_vcf
                ]
                
                # Add --full-annotation flag if requested
                if full_vcf_annotation:
                    vcf_args.append('--full-annotation')
                    # If we have separate features file, use it instead
                    if output_features and 'features_txt_output' in outputs:
                        # We need to merge predictions with features for full annotation
                        # For now, just note this in the log
                        self.logger.info("Full VCF annotation requested: all features and SHAP values will be added to VCF")
                    else:
                        self.logger.info("Full VCF annotation requested: all features and SHAP values will be added to VCF")
                
                # Add original command if available
                if self.original_command:
                    vcf_args.extend(['--command', self.original_command])
                
                self._run_step(
                    'add_NMDannot_to_vcf.py',
                    vcf_args,
                    "8. Add NMD annotation to VCF"
                )
                outputs['final_vcf'] = str(final_vcf)
        
        self.logger.info("=" * 70)
        self.logger.info("predNMD Pipeline Completed Successfully!")
        self.logger.info("=" * 70)
        
        if 'final_vcf' in outputs:
            self.logger.info(f"Final annotated VCF: {outputs['final_vcf']}")
        if 'predictions_txt' in outputs:
            self.logger.info(f"Predictions table: {outputs['predictions_txt']}")
        if 'features_txt_output' in outputs:
            self.logger.info(f"Features table: {outputs['features_txt_output']}")
        
        # Clean up intermediate files if requested
        if not keep_intermediate:
            self.logger.info("Cleaning up intermediate files...")
            intermediate_files = []
            
            # Define which files to clean based on keep_files
            if 1 not in keep_files and 'protein_coding_vcf' in outputs:
                intermediate_files.append(outputs['protein_coding_vcf'])
            
            # VEP files cleanup logic
            if 2 not in keep_files:
                # Clean filtered VCF if we went past step 3
                if 'vep_filtered_vcf' in outputs and end_step > 3:
                    intermediate_files.append(outputs['vep_filtered_vcf'])
                # Clean complete VEP annotated VCF only if step 8 is complete or not run
                # (we need it for step 8, so only clean after step 8 completes)
                if 'vep_annotated_vcf' in outputs and end_step > 8:
                    intermediate_files.append(outputs['vep_annotated_vcf'])
            
            # Always clean intermediate feature files from steps 3 and 6 (unless --keep-all)
            if 'features_txt' in outputs and end_step > 6:
                intermediate_files.append(outputs['features_txt'])
            if 'translationai_fasta' in outputs and end_step > 6:
                intermediate_files.append(outputs['translationai_fasta'])
            if 'features_translationai_txt' in outputs and end_step > 7:
                intermediate_files.append(outputs['features_translationai_txt'])
            
            # Always clean up VEP and TranslationAI intermediate files
            vep_summary = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz.vep.vcf_summary.html"
            vep_warnings = self.output_dir / f"{self.sample_name}.protein_coding.vcf.gz.vep.vcf_warnings.txt"
            translationai_h5 = self.output_dir / f"{self.sample_name}.translationAI.h5"
            translationai_orfs = self.output_dir / f"{self.sample_name}.translationAI_predORFs_0.0001_0.0001.txt"
            translationai_tis = self.output_dir / f"{self.sample_name}.translationAI_predTIS_0.0001.txt"
            translationai_tts = self.output_dir / f"{self.sample_name}.translationAI_predTTS_0.0001.txt"
            
            for f in [vep_summary, vep_warnings, translationai_h5, translationai_orfs, 
                     translationai_tis, translationai_tts]:
                if os.path.exists(f):
                    intermediate_files.append(str(f))

            for f in intermediate_files:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                        self.logger.debug(f"Removed: {f}")
                    except OSError as e:
                        self.logger.warning(f"Could not remove {f}: {e}")
        
        return outputs
    
    def run_step(self, step_number: int, input_file: str, **kwargs) -> str:
        """
        Run a single step of the pipeline
        
        Args:
            step_number: Step number (1-8)
            input_file: Input file for this step
            **kwargs: Additional arguments for the step
        
        Returns:
            Output file path
        """
        step_methods = {
            1: self._step1_filter_protein_coding,
            2: self._step2_vep_annotation,
            3: self._step3_add_features,
            4: self._step4_add_loeuf,
            5: self._step5_add_phylop,
            6: self._step6_add_translationai,
            7: self._step7_rf_prediction,
            8: self._step8_add_to_vcf
        }
        
        if step_number not in step_methods:
            raise ValueError(f"Invalid step number: {step_number}. Must be 1-8.")
        
        return step_methods[step_number](input_file, **kwargs)
    

