"""
scQPAS: Single-cell Quantification of PolyAdenylation Sites

A bioinformatics tool for analyzing polyadenylation sites in single-cell RNA-seq data.
"""

__version__ = "0.1.0"
__author__ = "Zavolan Lab"

# Import main API for convenient access
from .cli import main
from .core import run_pipeline

__all__ = [
    "main",
    "run_pipeline",
    "__version__",
    "__author__",
]
