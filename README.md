# SNPDash

## Overview

SNPDash is a web-based tool designed to semi-automate the process of prioritizing candidate genes associated with trait-linked SNPs. This interactive app helps to efficiently identify and annotate potential causal genes by combining GWAS results with gene functional annotation data.

## Features

- Manhattan plot for trait-SNP associations
- Genome-wide significance thresholds
- Click-to-select SNP functionality
- Status tracking (Todo, Done, Skipped)

## Installation

1. Create a new conda environment:
```bash
conda create -n snp-dash python=3.9
conda activate snp-dash
```

2. Clone the repository:
```bash
git clone [your-repository-url]
cd SNP-dash
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage

1. Start the application:
```bash
python app.py
```

2. Open your web browser and navigate to:
```
http://127.0.0.1:8050/
```

## Data Management

### Input Data
The application requires the following input files in the `data` directory:
- `snp_trait_data.csv`: GWAS results with trait-SNP associations
- `gene_data.csv`: Gene annotations and functional information

### Output Data
The application maintains two types of output:

1. **Pickle File** (`annotations.pkl`):
   - Automatically saves progress

2. **CSV Export**:
   - Generates a comprehensive report of prioritized genes