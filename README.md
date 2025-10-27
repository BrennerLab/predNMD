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

NMD requires several annotation files. Download and configure paths in `config.yaml`:

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

NMD runs the following steps if starting from an unannotated VCF file:

1. **Filter VCF**: Keep only variants located in protein-coding regions
2. **VEP Annotation**: Annotate variants with VEP 
3. **Add PTC Features**: Add features for each variant, which will be input to the Random Forest model
4. **Add LOEUF/PhyloP**: Add LOEUF and phyloP scores to the feature table
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
- `dis_to_first_outframeAUG`: distance of PTC to the first downstream out-of-frame AUG
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

## All command options
```bash
deepnmd run -h
usage: deepnmd run [-h] -i INPUT -o OUTPUT_DIR -s SAMPLE_NAME [-c CONFIG] [--skip-filtering] [--gene GENE]
                   [--skip-vep] [--from-features] [--no-vcf-output] [--only-vcf-annotation] [--predictions-file FILE]
                   [--keep-all] [--keep-filtered-vcf] [--keep-features] [--gtf GTF] [--genome GENOME] [--cds CDS]
                   [--gnomad GNOMAD] [--phylop PHYLOP] [--m6a M6A] [--expression EXPRESSION] [--vep-path VEP_PATH]
                   [--vep-cache VEP_CACHE] [--assembly {GRCh37,GRCh38}] [--model-dir MODEL_DIR] [--threads THREADS]
                   [--no-canonical]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file (VCF, VEP-annotated VCF, or feature table)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory
  -s SAMPLE_NAME, --sample-name SAMPLE_NAME
                        Sample name
  -c CONFIG, --config CONFIG
                        Configuration YAML file

Workflow options:
  --skip-filtering      Skip protein-coding region filtering (input already filtered)
  --gene GENE           Filter variants to specific gene (by SYMBOL or Ensembl ID) - skips Step 1
  --skip-vep            Skip VCF filtering and VEP annotation (input already VEP-annotated)
  --from-features       Start from existing feature table (input is features.txt)
  --no-vcf-output       Generate predictions table only (skip VCF annotation)
  --only-vcf-annotation
                        Only add predictions to VCF (requires --predictions-file)
  --predictions-file FILE
                        Pre-computed predictions file (for --only-vcf-annotation)

File retention options:
  --keep-all            Keep all intermediate files
  --keep-filtered-vcf   Keep protein-coding filtered VCF
  --keep-features       Keep feature extraction outputs

Reference files (override config):
  --gtf GTF             GTF annotation file
  --genome GENOME       Genome FASTA file
  --cds CDS             CDS FASTA file

Annotation files (override config):
  --gnomad GNOMAD       gnomAD constraint metrics file
  --phylop PHYLOP       phyloP bigWig file
  --m6a M6A             m6A annotation file
  --expression EXPRESSION
                        Gene expression file

VEP options (override config):
  --vep-path VEP_PATH   Path to VEP executable
  --vep-cache VEP_CACHE
                        VEP cache directory
  --assembly {GRCh37,GRCh38}
                        Genome assembly

Model options (override config):
  --model-dir MODEL_DIR
                        Directory containing Random Forest models

Runtime options:
  --threads THREADS     Number of threads
  --no-canonical        Include non-canonical transcripts
```


## Interpretation Guide

### Mechanism Classification

**N_terminal (N-terminal rescue):**
- Variant creates stop codon early in transcript
- Downstream in-frame AUG allows reinitiation
- Produces truncated protein missing N-terminus

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
  

## Citation

If you use NMD in your research, please cite:

```
[Paper citation here TODO]
```

## License

MIT License - see LICENSE file for details


## Contact

- Issues: [GitHub Issues](https://github.com/yaqisu/NMD/issues)
- Email: yaqisu@berkeley.edu

