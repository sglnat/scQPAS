# scQPAS - Single-Cell Quantification of PolyAdenylation Sites

A Python package for detecting and quantifying distances between reads and polyadenylation sites in single-cell RNA-seq data.

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

# Install scQPAS
pip install .
```

### Development Installation

For development work with editable mode:

```bash
git clone https://github.com/sglnat/scQPAS.git
cd scQPAS

# Create conda environment with all dependencies
conda env create -f environment.yml

# Activate the environment
conda activate scqpas

# Install scQPAS in editable mode
pip install -e .
```

### Command-line usage

Once installed and the conda environment is activated, use the CLI:

```bash
scqpas --bam sample.bam --gtf annotation.gtf --output results.csv
```

The polyadenylation site (PAS) is now **automatically detected** from the read data by identifying
the cleavage position with the most supporting evidence. No manual PAS specification is needed.

For all options:
```bash
scqpas --help
```

### Programmatic usage

```python
from scqpas import run_pipeline

run_pipeline(
    bam_path='sample.bam',
    gtf_path='annotation.gtf',
    output_path='results.csv'
    # Note: PAS is automatically detected from reads
)
```

## Configuration

scQPAS uses YAML configuration files for managing pipeline parameters. You can use the default configuration or provide your own.

### Using Default Configuration

The default configuration is located at `scqpas/config/defaults.yaml` and includes:
- **PolyA detection**: percentage threshold (80%), length threshold (5 bp)
- **Output settings**: default filenames and BEDtools score
- **Logging**: default level (INFO) and log file location

### Using Custom Configuration

Create a custom YAML config file with your desired parameters:

```yaml
# my_config.yaml
polya_detection:
  percentage_threshold: 85
  length_threshold: 6
  short_polyA_length_cutoff: 5
  short_polyA_required_percentage: 100
  use_fixed_coordinates: false

output:
  default_output_file: "my_results.csv"
  bedtools_score: 1000
  column_separator: "\t"

logging:
  default_level: "DEBUG"
  default_file: "my_pipeline.log"
  enable_file_logging: true
```

Then use it with the CLI:

```bash
scqpas --config my_config.yaml --bam sample.bam --gtf annotation.gtf --pas pas_sites.bed
```

### Parameter Priority

Command-line arguments override config file values, which override defaults:

```
Default config < Custom config file < CLI arguments
```

For example:
```bash
# Uses config file value for percentage_threshold
scqpas --config my_config.yaml --bam sample.bam --gtf annotation.gtf --pas pas_sites.bed

# CLI argument overrides config file
scqpas --config my_config.yaml --bam sample.bam --gtf annotation.gtf --pas pas_sites.bed --percentage-threshold 90
```

## Project Status

**Current Stage:**

- ✅ Modularized Python package with proper packaging structure
- ✅ CLI interface via Click
- ✅ Core pipeline functionality for single CPA site analysis
- ✅ BAM extraction with polyA detection
- ✅ GTF parsing and feature extraction
- ✅ CIGAR-based intron validation
- ✅ Distance calculation from CPA sites
- ✅ Logging framework implemented
- ✅ YAML-based configuration management
- ✅ Comprehensive docstrings (NumPy format) for all functions
- ✅ Full type hints on all function signatures

**To be implemented:**

- [ ] Multi-CPA site analysis
- [ ] Visualization of distance distributions
- [ ] Likelihood calculation methods
- [ ] Comprehensive unit test coverage


## Repository Structure

```
scQPAS/
├── README.md                          # This file
├── LICENSE                            # Apache 2.0 (chosen arbitrarly for now)
├── pyproject.toml                     # Package configuration
├── environment.yml                    # Conda environment specification
│
├── scqpas/                            # Main Python package
│   ├── __init__.py                    # Package initialization, version info
│   ├── cli.py                         # CLI interface (Click-based)
│   ├── core.py                        # Pipeline orchestration logic
│   ├── logging_config.py              # Logging configuration
│   ├── config_manager.py              # Configuration management
│   ├── extract_reads_BAM.py           # BAM file parsing, polyA detection
│   ├── extract_annotation_GTF.py      # GTF parsing, feature extraction
│   ├── extend_1kb_GTF.py              # GTF terminal exon extension
│   ├── extract_cigar_BAM.py           # CIGAR string parsing and validation
│   ├── calculate_distances.py         # Distance computation logic
│   ├── bedtools_intersections.py      # BEDtools intersection wrapper
│   │
│   └── config/                        # Configuration files
│       └── defaults.yaml              # Default pipeline parameters
│
├── tests/                             # Unit tests
│   └── __init__.py
│
└── .gitignore                         # Git ignore rules
```


## Pipeline Overview

The pipeline processes single-cell BAM files to quantify distances between read cleavage sites and polyadenylation sites:

1. **Read Extraction** (`extract_reads_BAM.py`): Extract reads from BAM, detect polyA tails by nucleotide composition
2. **Annotation Extraction** (`extract_annotation_GTF.py`): Parse GTF to generate genes, exons, and introns
3. **CIGAR Extraction** (`extract_cigar_BAM.py`): Parse CIGAR strings to derive intron coordinates from reads
4. **Bedtools Intersection** (`bedtools_intersections.py`): Intersect reads with annotated features
5. **CIGAR Filtering** (`extract_cigar_BAM.py`): Validate reads by comparing CIGAR-derived introns with annotation
6. **Distance Calculation** (`calculate_distances.py`): Compute distances from CPA sites to read ends

Additional:
**GTF Extension** (`extend_1kb_GTF.py`): Extend terminal exons by 1kb to capture 3' UTR regions. This file will be used for multi-PAS distance

### Data Flow

![Pipeline Diagram](dataflow_1.jpg)


## Requirements

- Python ≥ 3.10
- **Core dependencies**: pandas ≥ 1.5.0, numpy ≥ 1.23.0, click ≥ 8.0, pysam ≥ 0.21.0, pyyaml ≥ 6.0
- **External**: bedtools (for sequence intersection operations)

**Development dependencies** (optional, for `pip install -e ".[dev]"`):
- pytest ≥ 7.0 (testing)
- black ≥ 23.0 (code formatting)
- mypy ≥ 1.0 (static type checking)
- sphinx ≥ 5.0 (documentation generation)

See `pyproject.toml` for complete version specifications.


## Development

### Setup for development

```bash
git clone https://github.com/sglnat/scQPAS.git
cd scQPAS
pip install -e ".[dev]"
```

### Code quality

The codebase follows best practices for production-quality Python:

**Docstrings**: All functions have comprehensive NumPy-style docstrings describing:
- Purpose and behavior
- Parameters with types and descriptions
- Return values with types
- Examples where applicable

**Type Hints**: All function signatures include full type annotations:
```python
from typing import Optional, Union, Tuple
import pandas as pd
from .config_manager import ConfigManager

def extract_exons(
    gtf_file_path: str,
    exons_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
```

These enable:
- IDE autocomplete and inline documentation
- Static type checking with mypy: `mypy scqpas/`
- Better code clarity and maintainability
- Programmatic API exploration

**Code Formatting**: Use Black for consistent code style:
```bash
black scqpas/
```

### Code organization

- **`scqpas/`**: Core Python modules with importable functions
- **`scqpas/config/`**: Configuration files and config manager
- **`cli.py`**: CLI interface that wraps core functions
- **`tests/`**: Unit tests (to be expanded)
- Centralized config management via YAML files (see Configuration section above)


## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.