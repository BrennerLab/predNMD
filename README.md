# predNMD: prediction of nonsense-mediated mRNA decay (NMD) 

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**predNMD** is a comprehensive Python package for predicting if given stop-gain variants would trigger NMD. It uses a Random Forest model trained from various data sources, and for those variants predicted to not trigger NMD, the software can further predict if the variant would lead to N- or C-terminal truncated protein utilizing SHAP analysis.


## Installation

### Install predNMD

```bash
#It’s highly recommended to start a fresh environment to avoid potential dependency conflicts
conda create -n predNMD python=3.10 
conda activate predNMD

git clone https://github.com/BrennerLab/predNMD.git
cd predNMD
pip install .
```

### Download required Data Files

predNMD requires several annotation files. Download and configure paths in `config.yaml`:

1. LOEUF scores (can be downloaded with download_data.py)
2. PhyloP conservation scores (can be downloaded with download_data.py)
3. m6A modification data (included in /data)
4. Gene expression data (included in /data)
5. Ensembl reference genome (FASTA), reference GTF, and reference CDS (FASTA) (can be downloaded with download_data.py)
6. Pre-trained Random Forest model (included in /model)

```bash
# To see what datasets are available for download:
python download_data.py --list

# e.g., Download the LOEUF scores, hg38 phyloP conservation scores, and all three GRCh38 Ensembl reference files (release 104)
python download_data.py --datasets gnomad phylop-hg38 ensembl-all --assembly GRCh38 --ensembl-release 104

# You can also use download_data.py to download the Ensembl VEP cache, e.g.,
python download_data.py --datasets ensembl-vep --assembly GRCh38 --ensembl-release 104 --data-dir {OUTPUT_DATA_DIR}
tar -xzf {OUTPUT_DATA_DIR}/homo_sapiens_vep_104_GRCh37.tar.gz # decompress the downloaded cache
```

### Install TranslationAI 

**NOTE:** TranslationAI should be installed after predNMD has been installed.

```bash
git clone https://github.com/rnasys/TranslationAI.git
cd TranslationAI
python setup.py install
pip install tensorflow #TranslationAI dependency
```

### (Optional) The following stand-alone tools are only required if your input is raw VCF without Ensembl VEP annotation, while unnecessary if your input is Ensembl VEP annotated file: 

1. bedtools
```bash
conda install bioconda::bedtools
```

2. bcftools
```bash
conda install bioconda::bcftools
```

3. Ensembl VEP: please refer to [Ensembl VEP documentation](http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html) for guidance of downloading and installing Ensembl VEP


## Configuration

Edit `config.yaml` to set paths (need to be **absolute path**, not relative path) to required data files:

```yaml
# Ensembl reference files (REQUIRED)
reference:
  gtf_file: /path/to/reference.gtf # can be downloaded via download_data.py
  genome_fasta: /path/to/genome.fa # can be downloaded via download_data.py
  cds_fasta: /path/to/cds.fa # can be downloaded via download_data.py

# Annotation databases (REQUIRED)
annotation:
  gnomad_file: /path/to/gnomad.constraint_metrics.tsv # can be downloaded via download_data.py
  phylop_bigwig: /path/to/phyloP.bw # can be downloaded via download_data.py
  m6a_file: /path/to/m6A_annotations.txt # provided at /data/hg19_m6A-Atlas_highRes_all.txt.gz or /data/hg38_m6A-Atlas_highRes_all.txt.gz
  expression_file: /path/to/gene_expression.csv # provided at /data/GTEx_mean_expression_per_gene.csv

# Ensembl VEP configuration (REQUIRED if using Ensembl VEP)
vep:
  vep_path: /path/to/vep #path to vep executable
  cache_dir: /path/to/.vep #path to vep cache, can be downloaded via download_data.py
  assembly: GRCh37  # or GRCh38

# Machine learning model (REQUIRED)
model:
  model_dir: /path/to/RF_models/ # provided at /model

# Runtime settings
runtime:
  threads: 32 # change to the number of threads you want to use
  canonical_only: true # when set to "true", each variant will only be assigned to one transcript
                       # (with canonical transcript being prioritized), set to "false" if you want
                       # to include all potential isoforms containing the variants, which means
                       # one variant could correspond to multiple transcripts.
```


## Quick Start

- Basic usage
  ```bash
  prednmd run -i input.vcf -o output_dir -s output_prefix -c /path/to/config.yaml
  ```

- If your input is VEP-annotated VCF, please make sure it includes annotation for allele frequency. By default, predNMD will look for "gnomAD_AF" or "gnomADg_AF", but you could specify which allele frequency in your VEP-annotated input you want to use by adding the option `--af-column`, e.g., if you want to use African population allele frequency from gnomAD:
  ```bash
  prednmd run -i input.vcf -o output_dir -s output_prefix -c /path/to/config.yaml --af-column gnomADg_AFR_AF
  ```

- By default, predNMD will assign one transcript to each variant prioritized with canonical transcript, which means if your input is VEP-annotated VCF, it should contain the annotation for canonical transcripts (i.e., VEP should be run with `--canonical`). However, if you prefer to include all the transcript isoforms for each variant, you can run predNMD with `--no-canonical` (and in this case your VEP-annotated input file is not required to include the annotation for canonical transcripts):
  ```bash
  prednmd run -i input.vcf -o output_dir -s output_prefix -c /path/to/config.yaml --no-canonical
  ```
  
- If you only want to analyze variants in a specific gene, you can specify the gene symbol or the Ensembl gene ID with `--gene`, e.g., for variants in BRCA1:
  ```bash
  prednmd run -i input.vcf -o output_dir -s output_prefix -c /path/to/config.yaml --gene BRCA1 # or --gene ENSG00000012048
  ```
  
- By default, predNMD checks whether a SNV truly introduces a premature stop codon, even if it is annotated as “stop_gained” by Ensembl VEP, and will also check whether the reference allele matches at the variant position. This helps prevent mismatches caused by differences between the reference genome and GTF used by VEP and those provided to predNMD. However, you can use `--skip-ptc-check' to skip these verifications:
  ```bash
  prednmd run -i input.vcf -o output_dir -s output_prefix -c /path/to/config.yaml --skip-ptc-check
  ```

- If you also want to output the feature table containing the calculated features for input variants in addition to the prediction results, you can use the `--output-features` option:
  ```
  prednmd run -i input.vcf -o output_dir -c /path/to/config.yaml --output-features
  ```

- By default, predNMD will only append the prediction results to the INFO field of the VCF file. If you also want to annotate the VCF file with all the calculated features as well as SHAP values, you can use the `--full-vcf-annotation` option:
  ```
  prednmd run -i input.vcf -o output_dir -c /path/to/config.yaml --full-vcf-annotation
  ```

- If you already have pre-computed features that can directly serve as input to the RF model (please refer to [Input Requirements](https://github.com/yaqisu/NMD/tree/main?tab=readme-ov-file#input-requirements) below for the format of the feature table), you can run predNMD with `--from-features`
  ```bash
  prednmd run -i features.txt -o output_dir -c /path/to/config.yaml --from-features
  ```

- If you want to initialize a template config file
  ```bash
  prednmd init-config -o config.yaml
  ```


## Alternative option: using pre-built Docker image 

If any of the installation step fails and cannot be resolved, you can also directly use our pre-built Docker image as follows:

1. Pull the Docker image
```bash
docker pull brennerlab/prednmd:v1.0.0
```

2. Gather all your data files (including input and all reference files) into one directory, assuming it's {DATA_DIR} for the following example

3. Under {DATA_DIR} that contains all your data files, create a `config.yaml` as below. Please replace the {} parts with the actual **file name** under your {DATA_DIR}.
```yaml
reference:
  gtf_file: /data/{Homo_sapiens.GRCh37.87.gtf}
  genome_fasta: /data/{Homo_sapiens.GRCh37.dna.primary_assembly.fa}
  cds_fasta: /data/{GRCh37.CDS.fa}

annotation:
  gnomad_file: /data/{gnomad.v4.1.constraint_metrics.tsv}
  phylop_bigwig: /data/{hg19.100way.phyloP100way.bw}
  m6a_file: /app/data/hg19_m6A-Atlas_highRes_all.txt.gz # or change to /app/data/hg38_m6A-Atlas_highRes_all.txt.gz if using hg38 assembly
  expression_file: /app/data/GTEx_mean_expression_per_gene.csv

vep:
  vep_path: /opt/ensembl-vep-release-104.3/vep 
  cache_dir: /VEP_cache 
  assembly: GRCh37  # or change to GRCh38

model:
  model_dir: /app/model

runtime:
  threads: 32 # change to the number of threads you want to use on your local machine
  canonical_only: true # set to false if you want to include all transcripts
```

4. Under {DATA_DIR}, run prednmd via docker, e.g.,:
```bash
# If your input is VCF file unannotated by Ensembl VEP
docker run \
-v $(pwd):/data \ # Mount your local data path to the data path within the docker container
-v {PATH_TO_VEP_CACHE}:/VEP_cache \ # Replace {PATH_TO_VEP_CACHE} with the absolute path to VEP cache downloaded on your local machine
brennerlab/prednmd:v1.0.0 \
prednmd run -i /data/{YOUR_INPUT_VCF} -s {OUTPUT_PREFIX} -o /data/{OUTPUT_DIR_NAME} -c /data/config.yaml # you can specify different options as you need

# If your input is VEP-annotated VCF
docker run \
-v $(pwd):/data \ # Mount your local data path to the data path within the docker container
brennerlab/prednmd:v1.0.0 \
prednmd run -i /data/{YOUR_INPUT_VCF} -s {OUTPUT_PREFIX} -o /data/{OUTPUT_DIR_NAME} -c /data/config.yaml # you can specify different options as you need
```


## Output Files

### Predictions File (`{OUTPUT_PREFIX}.with_predictions.txt`)

Tab-delimited file with the following columns:

**Variant Information:**
- `CHR`: Chromosome containing the variant
- `POS`: Position of the variant on the chromsome
- `REF_ALLELE`: Reference allele at the variant position
- `ALT_ALLELE`: Alternative allele, i.e., the variant itself
- `transcript_id`: Ensembl ID of transcript containing the variant
- `gene_id`: Ensembl ID of gene containing the variant

**Core Prediction:**
- `nmd_trigger_probability`: Probability of triggering NMD (0-1)

**Mechanism Classification** (for NMD-not trigger cases only, i.e., when nmd_trigger_probability < 0.5):
- `mechanism_classification`: N_terminal, C_terminal, or Uncertain
- `n_terminal_probability`: Probability of N-terminal rescue mechanism
- `c_terminal_probability`: Probability of C-terminal rescue mechanism

### Annotated VCF (`{OUTPUT_PREFIX}.NMDannot.vcf`)

VCF file with added INFO fields:
- `NMD_PROB`: NMD-trigger probability
- `MECH_CLASS`:  Mechanism classification for variants predicted as not triggering NMD (N_terminal, C_terminal, or Uncertain)
- `N_TERMINAL_PROB`: N-terminal truncation probability
- `C_TERMINAL_PROB`: C-terminal truncation probability

### Log files

- `{OUTPUT_PREFIX}_prednmd.log`: standard output during running predNMD
- `{OUTPUT_PREFIX}.feature_added.warnings.log`: all warning messages during feature calculation for each variant (e.g., potential reference allele mismatch, etc.)

## Input Requirements

### VCF Input

- Can be VCF either unannotated (will run Ensembl VEP automatically) or annotated with Ensembl VEP (must include allele frequency annotation)
- All input alleles must be provided in the **forward-strand orientation**, regardless of whether the PTC generated by the allele occurs on a forward- or reverse-stranded transcript.

### Feature Table Input

If you would like to directly start from model prediction step, please provide a tab-delimited file with the following features:

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


## Pipeline Steps

NMD runs the following steps if starting from an unannotated VCF file:

1. **Filter VCF**: Keep only variants located in protein-coding regions
2. **Ensembl VEP Annotation**: Annotate variants with Ensembl VEP 
3. **Add PTC Features**: Add features for each variant, which will be input to the Random Forest model
4. **Add LOEUF/PhyloP**: Add LOEUF and phyloP scores to the feature table
5. **TranslationAI**: Apply TranslationAI to get predicted TIS/TTS scores for downstream inframe AUG/PTC
6. **RF Prediction**: Apply Random Forest model with SHAP analysis 
7. **Annotate VCF**: Add prediction results to the INFO field of VEP-annotated VCF 


## All command options
```bash
prednmd run -h
usage: prednmd run [-h] -i INPUT -o OUTPUT_DIR -s SAMPLE_NAME [-c CONFIG] [--skip-filtering] [--gene GENE]
                   [--from-features] [--no-vcf-output] [--only-vcf-annotation] [--predictions-file FILE]
                   [--output-features] [--full-vcf-annotation] [--skip-ptc-check] [--af-column COLUMN]
                   [--keep-all] [--keep-filtered-vcf] [--gtf GTF_FILE] [--genome GENOME_FILE] [--cds CDS_FILE]
                   [--gnomad GNOMAD_FILE] [--phylop PHYLOP_FILE] [--m6a M6A_FILE]
                   [--expression EXPRESSION_FILE] [--vep-path VEP_PATH] [--vep-cache VEP_CACHE]
                   [--assembly {GRCh37,GRCh38}] [--model-dir MODEL_DIR] [--threads THREADS] [--no-canonical]

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
  --skip-filtering      Skip protein-coding region filtering
  --gene GENE           Filter variants to specific gene (by SYMBOL or Ensembl ID)
  --from-features       Start from existing feature table (tab-delimited)
  --no-vcf-output       Generate predictions table only (skip adding NMD results to INFO field of the VCF)
  --only-vcf-annotation
                        Just add NMD predictions to VCF (requires --predictions-file)
  --predictions-file FILE
                        Pre-computed predictions file (for --only-vcf-annotation)

Output options:
  --output-features     Output a separate feature table with all features and SHAP values 
  --full-vcf-annotation
                        Add all features and SHAP values to VCF INFO field in addition to standard NMD annotations

Feature extraction options:
  --skip-ptc-check      Skip PTC check for SNVs and also skip check for reference allele matching in step 3. 
  --af-column COLUMN    Specify which AF column to use in step 3 (default: auto-detect gnomAD_AF or gnomADg_AF)

File retention options:
  --keep-all            Keep all intermediate files
  --keep-filtered-vcf   Keep protein-coding filtered VCF

Reference files (override config):
  --gtf GTF_FILE        GTF annotation file
  --genome GENOME_FILE  Genome FASTA file
  --cds CDS_FILE        CDS FASTA file

Annotation files (override config):
  --gnomad GNOMAD_FILE  gnomAD constraint metrics file
  --phylop PHYLOP_FILE  phyloP bigWig file
  --m6a M6A_FILE        m6A annotation file
  --expression EXPRESSION_FILE
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
  

## Citation

If you use predNMD in your research, please cite:

Su Y & Brenner SE. predNMD: prediction of nonsense-mediated mRNA decay for improved clinical variant pathogenicity classification. Preprint at bioRxiv. <https://doi.org/10.64898/2026.06.20.733449> (2026).





