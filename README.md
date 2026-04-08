# scQPAS - Single-Cell Quantification of PolyAdenylation Sites

A Python package for quantifying distances between reads and polyadenylation sites in single-cell RNA-seq data.

## Overview

scQPAS analyzes single-cell RNA-seq BAM files to identify where transcripts are being cleaved and polyadenylated. For each read, it calculates the distance between the actual read start and known or inferred polyadenylation sites (PAS), providing insights into 3' end processing and transcript cleavage patterns across cell types.

The tool divides reads from an **input BAM file** into two categories: reads with a polyA tail (polyA reads) and reads without one (non-polyA reads). polyA read detection is based on the workflow [SCINPAS](https://academic.oup.com/nar/article/53/D1/D197/7893321). scQPAS then works on read set (RS) level: A RS consists of all reads with the same cell barcode and unique molecular identifier, with polyA RS containing at least one polyA read and non-polyA RS containing only non-polyA reads. Since each RS is sequenced from the same RNA molecule, all reads within belong to the same PAS. scQPAS enables transcriptomic distance calculation of each read to its PAS, using an **input BED file from [PolyASite Atlas v3.0](https://polyasite.unibas.ch/atlas_sc)** to include all possible PAS for a RS and an **input GTF file** to account for different splicing isoforms that a RS can belong to.

This workflow results in two files:
- **candidate read-PAS distances**: a file containing all reads and the distances of each read to all its possible PAS, respective of all possible transcript isoforms for this RS
- **ground truth distances**: a file containing the distances of all reads within a polyA RS to the inferred cleavage and polyadenylation site(s) (CPA) of this RS

When comparing the histogram plots of the ground truth distances and the candidate distances, the candidate reads belonging to the overlapping part of the histograms can be assigned to the PAS that the CPA of the ground truth data belongs to. **This is to be implemented computationally via Maximum Likelihood Calculation.**

**IMPORTANT:** Due to high-complexity intermediate data generation, it is recommended to only apply scQPAS to regions of interest that span a gene.

## Installation

### Prerequisites

**Miniconda or Anaconda** must be installed. This provides:
- The conda package manager
- Python 3.10+
- Ability to manage system dependencies (bedtools, samtools) and Python packages

### User Installation

```bash
git clone https://github.com/sglnat/scQPAS.git
cd scQPAS

# Create conda environment with all dependencies
conda env create -f environment.yml

# Activate the environment
conda activate scqpas
```

Installation of scQPAS depending on regular or editable mode:
```bash
# Install scQPAS in regular mode
pip install .

# Install scQPAS for development work in editable mode:
pip install -e .
```


### Command-line usage

Once installed and the conda environment is activated, use the CLI:

```bash
scqpas --bam sample.bam --gtf annotation.gtf --pas atlas.bed --output results.csv
```

For ground truth data, the polyadenylation site (PAS) is **automatically detected** from the read data by identifying the cleavage position with the most supporting evidence. For candidate distances, manual PAS specification is needed.

To relieve computational burden, it is recommended to use PAS input file prefiltered to stringency and region of interest. Alternatively, it is also possible to filter by region and stringency directly in scQPAS, using dedicated flags.

**IMPORTANT:** It is strongly recommended to use `--region` flag.

For all options:
```bash
scqpas --help
```


## Configuration

scQPAS uses YAML configuration files for managing pipeline parameters. The default configuration is at `scqpas/config/defaults.yaml`.

Create a custom YAML file to override defaults:

```yaml
# config.yaml
polya_detection:
  percentage_threshold: 85
  length_threshold: 6
```

Then pass it to the CLI:

```bash
scqpas --config config.yaml --bam sample.bam --gtf annotation.gtf --pas atlas.bed --output results.csv
```


## How It Works

1. **Read Extraction**: Extract reads from BAM, detect polyA tails
2. **Annotation Parsing**: Extract genes, exons, and introns from GTF
3. **CIGAR Processing**: Parse CIGAR strings to derive intron coordinates
4. **Feature Intersection**: Intersect reads with annotated features using bedtools
5. **Read Validation**: Filter reads using CIGAR-derived introns
6. **Distance Calculation**: Compute distances from polyadenylation sites to read ends

### Data Flow

![Pipeline Diagram](dataflow_1.jpg)



### Code Quality

- **Type hints**: Full type annotations on all functions enable IDE support and static checking
- **Docstrings**: NumPy-style docstrings for all functions
- **Formatting**: Use Black: `black scqpas/`
- **Type checking**: Run mypy: `mypy scqpas/`

### Structure

- `scqpas/`: Core modules
- `scqpas/config/`: Configuration files
- `cli.py`: CLI interface



## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.