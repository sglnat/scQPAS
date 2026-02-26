"""
CLI interface for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module provides the command-line interface for the scQPAS pipeline.
Core pipeline logic is in core.py.
"""

import sys
import os
import logging

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

import click
from src.core import run_pipeline
from src.logging_config import configure_logging

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    '--bam',
    type=click.Path(exists=True),
    required=True,
    help='Path to input BAM file (must be indexed)'
)
@click.option(
    '--gtf',
    type=click.Path(exists=True),
    required=True,
    help='Path to GTF annotation file'
)
@click.option(
    '--chr',
    type=str,
    required=True,
    help='Target chromosome (e.g., chr12)'
)
@click.option(
    '--pas',
    type=int,
    required=True,
    help='Polyadenylation site position (genomic coordinate)'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='distances.csv',
    help='Output file path (CSV format). Default: distances.csv'
)
@click.option(
    '--percentage-threshold',
    type=int,
    default=80,
    help='Min %% of A nucleotides in polyA region (0-100). Default: 80'
)
@click.option(
    '--length-threshold',
    type=int,
    default=5,
    help='Min length of polyA tail (bp). Default: 5'
)
@click.option(
    '--use-fc',
    is_flag=True,
    default=False,
    help='Use fixed soft-clipped coordinates from BAM XO/XF tags'
)
@click.option(
    '--log-level',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR'], case_sensitive=False),
    default='INFO',
    help='Logging level. Default: INFO'
)
@click.option(
    '--log-file',
    type=click.Path(),
    default='scqpas.log',
    help='Log file path. Default: scqpas.log (set to empty string to disable file logging)'
)
def main(bam, gtf, chr, pas, output, percentage_threshold, length_threshold, use_fc, log_level, log_file):
    """
    Calculate distances from CPA sites to reads.
    
    Processes a BAM file against a GTF annotation in 6 steps:
    1. Extract polyA-containing reads from BAM
    2. Extract annotation (genes, exons, introns) from GTF
    3. Extract CIGAR-derived introns from reads
    4. Intersect reads with annotated features (bedtools)
    5. Filter reads by CIGAR string validation
    6. Calculate distances from cleavage to PAS
    
    All processing in-memory except temporary files for bedtools (auto-cleaned).
    
    Example:
    
        scqpas --bam sample.bam --gtf annotation.gtf --chr chr12 --pas 6538371 --output results.csv
    """
    
    try:
        # Configure logging based on user's options
        log_level_int = getattr(logging, log_level.upper())
        # If log_file is empty string, pass None to disable file logging
        log_file_arg = log_file if log_file else None
        configure_logging(log_level=log_level_int, log_file=log_file_arg)
        
        run_pipeline(
            bam_path=bam,
            gtf_path=gtf,
            chr=chr,
            pas=pas,
            output_path=output,
            percentage_threshold=percentage_threshold,
            length_threshold=length_threshold,
            use_fc=use_fc
        )
        
    except click.ClickException:
        raise
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        raise click.ClickException(str(e))
    except RuntimeError as e:
        logger.error(f"Pipeline error: {e}")
        raise click.ClickException(str(e))
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
