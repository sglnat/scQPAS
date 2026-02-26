# scQPAS - Single-Cell Quantification of PolyAdenylation Sites

A Python package for detecting and quantifying distances between reads and polyadenylation sites in single-cell RNA-seq data.

## Installation

### From source (development mode)

```bash
git clone https://github.com/sglnat/scQPAS.git
cd scQPAS
pip install -e .
```

### Command-line usage

Once installed, use the CLI:

```bash
scqpas --bam sample.bam --gtf annotation.gtf --chr chr12 --pas 6538371 --output results.csv
```

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
    chr='chr12',
    pas=6538371,
    output_path='results.csv'
)
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
│   ├── extract_reads_BAM.py           # BAM file parsing, polyA detection
│   ├── extract_annotation_GTF.py      # GTF parsing, feature extraction
│   ├── extend_1kb_GTF.py              # GTF terminal exon extension
│   ├── extract_cigar_BAM.py           # CIGAR string parsing and validation
│   ├── calculate_distances.py         # Distance computation logic
│   └── bedtools_intersections.py      # BEDtools intersection wrapper
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
- Dependencies: pandas, numpy, click, pysam, bedtools

See `pyproject.toml` for complete version specifications.


## Development

### Setup for development

```bash
git clone https://github.com/sglnat/scQPAS.git
cd scQPAS
pip install -e ".[dev]"
```

### Code organization

- **`scqpas/`**: Core Python modules with importable functions
- **`cli.py`**: CLI interface that wraps core functions
- **`tests/`**: Unit tests (to be expanded)
- No config file for parameters (all passed as CLI args)
- Nextflow files localization


## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.