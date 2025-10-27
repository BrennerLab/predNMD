# DeepNMD: Deep Learning for NMD Escape Prediction

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DeepNMD is a comprehensive Python package for predicting Nonsense-Mediated Decay (NMD) escape mechanisms from genetic variants. It uses Random Forest models with SHAP analysis to classify variants as N-terminal or C-terminal rescue mechanisms.

## Features

- **Complete Pipeline**: From VCF to predictions in one command
- **Gene Filtering**: Analyze variants in specific genes of interest
- **Flexible Input**: Start from VCF or pre-computed feature tables
- **SHAP Analysis**: Interpretable predictions with mechanism classification
- **VEP Integration**: Automated VEP annotation (or use pre-annotated VCF)
- **Multi-threaded**: Parallel processing for faster analysis

## Installation

### Prerequisites

- Python 3.8 or higher
- VEP (Variant Effect Predictor) installed and configured
- Reference genome (GRCh38)

### Install DeepNMD

```bash
# Install from PyPI (when published)
pip install deepnmd

# Or install from source
git clone https://github.com/yourusername/deepnmd.git
cd deepnmd
pip install -e .
```

### Required Data Files

DeepNMD requires several annotation files. Download and configure paths in `config.yaml`:

1. LOEUF scores
2. PhyloP conservation scores
3. m6A modification data
4. Gene expression data
5. Reference genome (FASTA)
6. Pre-trained Random Forest model

## Quick Start

### Basic Usage

```bash
# Run complete pipeline
deepnmd -v input.vcf -o output_dir -c config.yaml
```

### Gene-Specific Analysis

```bash
# Analyze variants in BRCA1 only
deepnmd -v input.vcf -o output_dir --gene BRCA1

# Using Ensembl Gene ID
deepnmd -v input.vcf -o output_dir --gene ENSG00000012048
```

### Advanced Usage

```bash
# Skip VEP annotation (already annotated)
deepnmd -v annotated.vcf -o output_dir --skip-vep

# Start from feature table (step 7: model prediction only)
deepnmd -f features.txt -o output_dir

# Get predictions without VCF output
deepnmd -v input.vcf -o output_dir --no-vcf-output

# Keep intermediate VEP file
deepnmd -v input.vcf -o output_dir --keep-vep-vcf

# Use multiple threads
deepnmd -v input.vcf -o output_dir -t 8
```

## Pipeline Steps

DeepNMD runs the following steps:

1. **Filter VCF**: Keep only stop_gained variants
2. **VEP Annotation**: Annotate variants with VEP (skippable)
3. **Extract PTC Features**: Extract premature termination codon features
4. **Add LOEUF/PhyloP**: Annotate with constraint and conservation scores
5. **TranslationAI**: Calculate translation efficiency scores
6. **Add m6A/Expression**: Annotate with modification and expression data
7. **RF Prediction**: Apply Random Forest model with SHAP analysis (starts here if using `--features`)
8. **Annotate VCF**: Add predictions back to VCF (optional)

## Output Files

### Predictions File (`nmd_predictions.txt`)

Tab-delimited file with the following columns:

**Core Predictions:**
- `nmd_trigger_probability`: Probability of triggering NMD (0-1)
- `nmd_escape_probability`: Probability of escaping NMD (0-1)

**Mechanism Classification** (for escape cases only):
- `mechanism_classification`: N_terminal, C_terminal, or Uncertain
- `n_terminal_probability`: Probability of N-terminal rescue mechanism
- `c_terminal_probability`: Probability of C-terminal rescue mechanism

**SHAP Contributions:**
- `n_terminal_escape_contrib`: Raw SHAP contribution from N-terminal features
- `c_terminal_escape_contrib`: Raw SHAP contribution from C-terminal features
- `general_escape_contrib`: Raw SHAP contribution from general features

**Relative Contributions:**
- `n_terminal_escape_relative`: Relative contribution (3-way split including general)
- `c_terminal_escape_relative`: Relative contribution (3-way split)
- `general_escape_relative`: Relative contribution (3-way split)
- `n_terminal_nt_ct_relative`: Relative contribution (N vs C only)
- `c_terminal_nt_ct_relative`: Relative contribution (N vs C only)
- `total_nt_ct_escape_score`: Combined N+C terminal escape score

### Annotated VCF (`nmd_annotated.vcf`)

VCF file with added INFO fields:
- `NMD_TRIGGER_PROB`: NMD trigger probability
- `NMD_ESCAPE_PROB`: NMD escape probability
- `NMD_MECHANISM`: Mechanism classification (N_terminal/C_terminal/Uncertain)
- `NMD_N_PROB`: N-terminal mechanism probability
- `NMD_C_PROB`: C-terminal mechanism probability

## Gene Filtering

The `--gene` option filters variants to a specific gene using VEP annotations:

### Gene Specification

You can specify genes by:
- **Symbol**: Human-readable gene name (e.g., `BRCA1`, `TP53`)
- **Ensembl ID**: Ensembl gene identifier (e.g., `ENSG00000012048`)

### How It Works

1. Filters VCF after stop_gained filtering but before further analysis
2. Matches against both VEP `SYMBOL` and `Gene` fields
3. Retains all transcripts for the specified gene
4. Fails with error if no variants found for the gene

### Examples

```bash
# Common cancer genes
deepnmd -v input.vcf -o brca1_analysis --gene BRCA1
deepnmd -v input.vcf -o tp53_analysis --gene TP53
deepnmd -v input.vcf -o pten_analysis --gene PTEN

# Using Ensembl IDs
deepnmd -v input.vcf -o analysis --gene ENSG00000141510  # TP53

# Gene filtering requires VCF input (not feature tables)
deepnmd -f features.txt -o output --gene BRCA1  # ERROR

# Combine with other options
deepnmd -v input.vcf -o output --gene BRCA1 --skip-vep -t 4
```

### Important Notes

- Gene filtering requires VEP-annotated VCF (or use with `--skip-vep` if already annotated)
- Cannot be used with `--features` input (feature tables don't contain gene information)
- Gene names are case-sensitive (use exact match)
- Returns error if gene not found in VCF

## Configuration

Edit `config.yaml` to set paths to required data files:

```yaml
# Reference data files
loeuf_file: /path/to/loeuf_scores.txt
phylop_file: /path/to/phylop_scores.txt
m6a_file: /path/to/m6a_annotations.txt
expression_file: /path/to/expression_data.txt
genome_file: /path/to/GRCh38.fa

# Model directory
model_dir: /path/to/rf_model

# Processing options
threads: 4
```

## Input Requirements

### VCF Input

- Must contain stop_gained variants
- Can be raw VCF (will run VEP) or VEP-annotated (use `--skip-vep`)
- GRCh38 assembly

### Feature Table Input

If starting from step 7 (model prediction), provide a tab-delimited file with these features:

**Categorical:**
- `50nt_rule`: Boolean
- `has_downstream_inframeAUG`: Boolean

**Continuous:**
- `CDS_position`, `dis_to_first_inframeAUG`, `dis_to_first_outframeAUG`
- `downstream_inframeAUG_translationAI`, `dis_to_exon_end`, `exon_length`
- `distance_to_stop`, `downstream_exons`, `dis_to_3utr_end`
- `CAI_25codon_upstream_diff`, `phyloP`, `upstream_exons`
- `AF`, `gc_content`, `LOEUF`, `PTC_translationAI`
- `Mean_Expression`, `m6A_CDS_length_normalized_unconstrained`
- `m6A_all_length_normalized_unconstrained`

## Examples

### Example 1: Complete Analysis

```bash
deepnmd -v variants.vcf -o results -c config.yaml -t 8
```

**Output:**
```
results/
├── filtered_stop_gained.vcf
├── vep_annotated.vcf (removed unless --keep-vep-vcf)
├── ptc_features.txt
├── features_with_loeuf_phylop.txt
├── features_with_tai.txt
├── complete_features.txt
├── nmd_predictions.txt
└── nmd_annotated.vcf
```

### Example 2: Gene-Specific Analysis

```bash
# Analyze BRCA1 variants only
deepnmd -v patient.vcf -o brca1_results --gene BRCA1 -t 4

# Multiple gene analyses
for gene in BRCA1 BRCA2 TP53 PTEN; do
    deepnmd -v patient.vcf -o ${gene}_analysis --gene $gene
done
```

### Example 3: Pre-annotated VCF

```bash
# Skip VEP if already annotated
deepnmd -v vep_annotated.vcf -o results --skip-vep
```

### Example 4: Feature Table Only

```bash
# Start from pre-computed features (step 7)
deepnmd -f my_features.txt -o predictions
```

### Example 5: Predictions Only (No VCF Output)

```bash
# Get prediction table without annotating VCF
deepnmd -v input.vcf -o results --no-vcf-output
```

## Interpretation Guide

### Mechanism Classification

**N_terminal (N-terminal rescue):**
- Variant creates stop codon early in transcript
- Downstream in-frame AUG allows reinitiation
- Produces truncated but potentially functional protein

**C_terminal (C-terminal rescue):**
- Variant follows 50-nucleotide rule (>50nt from last exon junction)
- Escapes NMD due to favorable exon structure
- Produces truncated protein missing C-terminus

**Uncertain:**
- Both mechanisms contribute equally
- Or neither mechanism shows clear positive contribution

### SHAP Contributions

SHAP values explain each prediction:
- **Positive values**: Push toward NMD escape
- **Negative values**: Push toward NMD trigger
- **Magnitude**: Importance of contribution

### Probability Interpretation

- `nmd_escape_probability > 0.5`: Predicted to escape NMD
- `n_terminal_probability > c_terminal_probability`: N-terminal mechanism more likely
- Higher probabilities indicate higher confidence

## Troubleshooting

### Common Issues

**1. Gene not found**
```
Error: No variants found for gene 'BRCA1'
```
- Check gene name spelling
- Verify gene has stop_gained variants in VCF
- Try Ensembl ID instead of symbol

**2. VEP annotation required**
```
Error: Gene filtering requires VEP annotation
```
- Run without `--skip-vep`, or
- Pre-annotate VCF with VEP manually

**3. Missing features**
```
Warning: Missing features: ['phyloP', 'LOEUF']
```
- Check config.yaml paths
- Ensure annotation files are in correct format

**4. Feature table with gene filter**
```
Error: Gene filtering requires VCF input
```
- Cannot use `--gene` with `--features`
- Filter VCF first, then extract features

### Getting Help

```bash
# Show help
deepnmd --help

# Enable verbose logging
deepnmd -v input.vcf -o output --verbose

# Save log to file
deepnmd -v input.vcf -o output --log-file analysis.log
```

## Citation

If you use DeepNMD in your research, please cite:

```
[Your paper citation here]
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Contact

- Issues: [GitHub Issues](https://github.com/yourusername/deepnmd/issues)
- Email: your.email@institution.edu

## Changelog

### Version 1.1.0 (2025-10-23)
- Added `--gene` option for gene-specific analysis
- Fixed step numbering (feature table now correctly starts at step 7)
- Improved error messages for gene filtering
- Added gene symbol and Ensembl ID support

### Version 1.0.0 (2025-10-15)
- Initial release
- Complete VCF-to-prediction pipeline
- SHAP-based mechanism classification
- Multi-threaded processing
