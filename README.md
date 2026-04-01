# Mutational Signatures Analysis

## Overview
This repository contains python script for analyzing mutational signatures in hpv-positive head and neck cancer samples. The pipeline uses SigProfilerMatrixGenerator to create mutation matrices from VCF files, SigProfilerExtractor to extract de-novo signatures, and SigProfilerAssignment to assign COSMIC signatures to the samples. 

The R script is used to analyze the sample structure based on the mutation signature structure using the output from SigProfilerAssignment (SBS96). We want to stratify the samples into subgroups based on the mutation signature structure and hopefully correlate with the cancer types (e.g., HPV-positive and HPV-negative). To-do: further analysis might involve statistical calculations to determine the significance of the subgroups, correlation with clinical data, identify the de novo signatures and decompose, what else(?). 


## Directory Structure
```
mutsig/
├── data/
│   ├── vcf/
│   │   ├── *.vcf
│   │   └── JOB_METADATA_SPM.txt
├── investigate_signatures.R
├── run_sigature_generate.py
├── requirements.txt
├── renv.lock
├── .Rprofile
├── renv/
│   └── activate.R
├── nv.lock
├── pyproject.toml
└── README.md
```

## Scripts

### `run_sigature_generate.py`
Automates the SigProfilerMatrixGenerator/SigProfilerExtractor/SigProfilerAssignment workflow:
1. Sanitizes VCF files (removes Windows line endings)
2. Installs the reference genome (GRCh37 by default)
3. Generates mutation matrices (SBS, DBS, ID)
4. Runs SigProfilerExtractor for de-novo signature extraction
5. Runs SigProfilerAssignment to assign COSMIC signatures to the samples.

**Usage:**
```bash
# Default (uses HN_HPV/HNhpv_vcf)
python run_sigature_generate.py

# Specify VCF directory
python run_sigature_generate.py --vcf-dir path/to/vcf

# Specify genome and project name
python run_sigature_generate.py --vcf-dir path/to/vcf --genome GRCh38 --project MyProject

# Skip signature extraction
python run_sigature_generate.py --vcf-dir path/to/vcf --skip-extraction
```

### `investigate_signatures.R`
Analyzes the sample structure based on the mutation signature structure using the output from SigProfilerAssignment. We want to stratify the samples into subgroups based on the mutation signature structure and hopefully correlate with the cancer types (HPV-positive and HPV-negative).

1. Loads `Assignment_Solution_Activities.txt`, the fitted cosmic 96 SBS signatures
2. Performs K-Means clustering on the samples based on the fitted signatures
3. Generates PCA plot (PC1 vs PC2)
4. Generates UMAP plot
5. Generates hierarchical clustering heatmap

**Usage:**
```bash
Rscript investigate_signatures.R
```

## Prerequisites
- Python 3.9
- requirements.txt for python packages
- R 4.3+ 
- renv.lock for R packages

## Installation

### Python Packages
```bash
pip install -r requirements.txt
```

### R Packages
```R
# Install renv if not already installed
if (!require("renv")) install.packages("renv")

# Initialize the project-local renv environment
renv::init()

# Install all packages listed in renv.lock
renv::restore()
```

## Example run

NOte: i have already tested the pipeline with the BRCA dataset and it works fine.

### Step 1: download the vcf files 

```bash
wget ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerAssignment/Example_data/BRCA.zip
```
This is the example dataset for the SigProfilerAssignment tool.

```bash
unzip BRCA.zip
```

This will create:
- BRCA/BRCA_vcf/ subdirectories where there are vcf files.

### Step 2: Run the SigProfiler pipeline 
```bash
python run_sigature_generate.py --vcf-dir BRCA/BRCA_vcf --min-sigs 2 --max-sigs 10 --nmf-reps 50
```

This will:
- generate the mutation matrices
- Run SigProfilerExtractor on the SBS96 matrix
- Extract de-novo signatures
- Map them to COSMIC reference signatures

### Step 3: Analyze Signature Structure
```bash
Rscript investigate_signatures.R <path/to/Assignment_Solution_Activities.txt>
```

The R script could be run as a script or interactively. It is better to run it with the interactive mode so to modify the number of clusters (k) if needed. 

This will generate:
- `PCA_Signatures.pdf` - PCA plot showing sample clustering
- `UMAP_Signatures.pdf` - UMAP plot showing sample clustering
- `Hierarchical_Clustering_Heatmap.pdf` - Heatmap of signature proportions

On BRCA dataset, the pipeline shows that the samples can be stratified into 3 or 4 subgroups based on the mutation signature structure.



