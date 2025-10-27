# DeepNMD: Deep Learning for NMD Escape Prediction

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**NMD(Name TBD)** is a comprehensive Python package for predicting if given nonsense or frameshift variants would trigger NMD. It uses a Random Forest model trained from various data sources, and for those variants predicted to not trigger NMD, the software can further predict if the variant would lead to N- or C-terminal truncated protein utilizing SHAP analysis.

## Installation

### Install NMD

```bash
#It’s highly recommended to start a fresh environment to avoid potential dependency conflicts
conda create -n NMD python=3.10 
conda activate NMD

git clone https://github.com/yaqisu/NMD.git
cd NMD
pip install .
```

### Download required Data Files

DeepNMD requires several annotation files. Download and configure paths in `config.yaml`:

1. LOEUF scores (can be downloaded with download_data.py)
2. PhyloP conservation scores (can be downloaded with download_data.py)
3. m6A modification data (already included in /data)
4. Gene expression data (already included in /data)
5. Ensembl reference genome (FASTA), reference GTF, and reference cDNA (FASTA)
6. Pre-trained Random Forest model (already included in /model)

```bash
# Download the gnomAD constraint metric (LOEUF scores) and phyloP conservation scores
python download_data.py

# You can also specify which dataset you want to download, e.g., 
python download_data.py --datasets gnomad phylop-hg38

# To see what datasets are available for download:
python download_data.py --list
```

### Install TranslationAI (should be installed after NMD has been installed)

```bash
git clone https://github.com/rnasys/TranslationAI.git
cd TranslationAI
python setup.py install
pip install tensorflow #TranslationAI dependency
```

### (Optional) The following stand-alone tools are only required if your input is raw VCF without VEP annotation, while unnecessary if your input is VEP annotated file: 

1. bedtools
```bash
conda install bioconda::bedtools
```

2. bcftools
```bash
conda install bioconda::bcftools
```

3. VEP: please refer to [VEP documentation](http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html) for guidance of downloading and installing VEP


## Quick Start

### Basic Usage

```bash
# If your input VCF has not been annotated by VEP
deepnmd run -i input.vcf -o output_dir -c /path/to/config.yaml
# If your input is VEP annotated VCF
deepnmd run -i vep_annotated.vcf -o output_dir -c /path/to/config.yaml --skip-vep
# If you want to initialize a template config file
deepnmd init-config -o config.yaml
```

### Gene-Specific Analysis

```bash
# For example, analyze variants in BRCA1 only (using gene symbol)
deepnmd run -i input.vcf -o output_dir -c /path/to/config.yaml --gene BRCA1

# Or you can also use Ensembl Gene ID
deepnmd run -i input.vcf -o output_dir -c /path/to/config.yaml --gene ENSG00000012048
```


## Pipeline Steps

DeepNMD runs the following steps:

1. **Filter VCF**: Keep only variants located in protein-coding regions
2. **VEP Annotation**: Annotate variants with VEP (skippable)
3. **Add PTC Features**: add features for each variant, which will be input to the Random Forest model
4. **Add LOEUF/PhyloP**: Annotate with constraint and conservation scores
5. **TranslationAI**: Apply TranslationAI to get predicted TIS/TTS scores for downstream inframe AUG/PTC
6. **RF Prediction**: Apply Random Forest model with SHAP analysis 
7. **Annotate VCF**: Add prediction results (probability of triggering NMD, probability of C-terminal truncation, probability of N-terminal truncation) back to the VCF 


## Output Files

### Predictions File (`{SAMPLE_NAME}.with_predictions.txt`)

Tab-delimited file with the following columns:

**Core Prediction:**
- `nmd_trigger_probability`: Probability of triggering NMD (0-1)

**Mechanism Classification** (for escape cases only):
- `mechanism_classification`: N_terminal, C_terminal, or Uncertain
- `n_terminal_probability`: Probability of N-terminal rescue mechanism
- `c_terminal_probability`: Probability of C-terminal rescue mechanism

**SHAP Contributions:**
- `n_terminal_escape_contrib`: Raw SHAP contribution from N-terminal features
- `c_terminal_escape_contrib`: Raw SHAP contribution from C-terminal features
- `general_escape_contrib`: Raw SHAP contribution from general features


### Annotated VCF (`{SAMPLE_NAME}.NMDannot.vcf`)

VCF file with added INFO fields:
- `NMD_PROB`: NMD trigger probability
- `N_TERMINAL_PROB`: N-terminal truncation probability
- `C_TERMINAL_PROB`: C-terminal truncation probability


## Configuration

Edit `config.yaml` to set paths to required data files:

```yaml
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
```

## Input Requirements

### VCF Input

- Can be either raw VCF (will run VEP) or VEP-annotated (use `--skip-vep`)

### Feature Table Input

If directly starting from model prediction step, please provide a tab-delimited file with these features:

**Categorical:**
- `50nt_rule`: True if the PTC is located <50nt upstream of the last exon-exon junction, False otherwise

**Continuous:**
- `CDS_position`: position of PTC on the CDS region
- `dis_to_first_inframeAUG`: distance of PTC to the first downstream inframe AUG
- `dis_to_first_outframeAUG: distance of PTC to the first downstream out-of-frame AUG
- `downstream_inframeAUG_translationAI`: translationAI score for the first downstream inframe AUG
- `dis_to_exon_end`: distance of PTC to the end of the PTC-containing exon
- `exon_length`: length of the PTC-containing exon
- `distance_to_stop`: distance of PTC to the normal stop codon
- `downstream_exons`: number of exons downstream of PTC
- `dis_to_3utr_end`: distance of PTC to the end of 3'UTR 
- `CAI_25codon_upstream_diff`: difference between CAI (Codon Adapation Index) 25 codons upstream of PTC and normal stop codon
- `phyloP`: phyloP score at the PTC site
- `upstream_exons`: number of exons upstream of PTC
- `AF`: allele frequency of the PTC-generating variant
- `gc_content`: GC content of the transcript 
- `LOEUF`: LOEUF score for the PTC-containing transcript
- `PTC_translationAI`: translationAI score for the PTC
- `Mean_Expression`: Gene mean expression taken from GTEx
- `m6A_CDS_length_normalized_unconstrained`: m6A counts between CDS start and PTC normalized by region length
- `m6A_all_length_normalized_unconstrained`: m6A counts across the whole transcript normalized by transcript length

## Examples

### Example 1: Complete Analysis

```bash
deepnmd -v variants.vcf -o results -c config.yaml 
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
