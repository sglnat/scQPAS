"""
CLI interface for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module provides the command-line interface for the scQPAS pipeline.
Core pipeline logic is in core.py.
"""

import logging
from typing import Optional

import click
from .core import run_pipeline
from .logging_config import configure_logging
from .config_manager import ConfigManager

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=True),
    default=None,
    help="Path to YAML config file (uses defaults if not provided)",
)
@click.option(
    "--bam",
    type=click.Path(exists=True),
    required=True,
    help="Path to input BAM file (must be indexed)",
)
@click.option(
    "--gtf",
    type=click.Path(exists=True),
    required=True,
    help="Path to GTF annotation file",
)
@click.option(
    "--pas",
    type=click.Path(exists=True),
    required=True,
    help="Path to polyadenylation sites BED file",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default=None,
    help="Output file path (CSV format). [default: distances.csv from config]",
)
@click.option(
    "--percentage-threshold",
    type=int,
    default=None,
    help="Min %% of A nucleotides in polyA region (0-100). [default: 80 from config]",
)
@click.option(
    "--length-threshold",
    type=int,
    default=None,
    help="Min length of soft-clipped region at 3' end (bp). [default: 5 from config]",
)
@click.option(
    "--terminal-exon-extension",
    type=int,
    default=None,
    help="Length in bp to extend terminal exons for PAS capture. [default: 1000 from config]",
)
@click.option(
    "--stringency",
    type=int,
    default=None,
    help="Minimum stringency level for PAS from sc PolyASite Atlas. [default: 80 from config]",
)
@click.option(
    "--region",
    type=str,
    default=None,
    help="Optional genomic region for PAS filtering (e.g., 'chr12:15671-26783' or 'chr12'). If not provided, all PAS are used.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default=None,
    help="Logging level [DEBUG|INFO|WARNING|ERROR]. Default: INFO (from config)",
)
@click.option(
    "--log-file",
    type=click.Path(),
    default=None,
    help="Log file path. Uses config default if not provided (set to empty string to disable)",
)
@click.option(
    "--debug-output",
    type=click.Path(),
    default=None,
    help="Directory to save intermediate pipeline files for debugging/inspection",
)
def main(
    config: Optional[str],
    bam: str,
    gtf: str,
    pas: str,
    output: Optional[str],
    percentage_threshold: Optional[int],
    length_threshold: Optional[int],
    terminal_exon_extension: Optional[int],
    stringency: Optional[int],
    region: Optional[str],
    log_level: Optional[str],
    log_file: Optional[str],
    debug_output: Optional[str],
) -> None:
    """
    Calculate distances from read starts to polyadenylation sites.

    Processes a BAM file against a GTF annotation and a PAS BED file to compute
    genomic distances from the 5' start of each read to the associated PAS position.

    The PAS site assignment is based on:
    1. Automatic detection of the most frequently used PAS in the dataset
    2. Association of PAS with transcripts they overlap
    3. Assignment of reads to transcripts and their associated PAS sites

    The distance output represents how far each read extends from its 5' start
    to the PAS position, with intron lengths subtracted for reads spanning introns.

    Note: The --pas argument specifies the BED file of known PAS regions but no
    longer forces a single fixed PAS position (which was the original behavior).

    All processing is in-memory except temporary files for bedtools (auto-cleaned).

    Configuration can be provided via:
    - Custom YAML config file (--config)
    - Defaults from scqpas/config/defaults.yaml
    - Command-line arguments (override config)

    Example:

        scqpas --bam sample.bam --gtf annotation.gtf --pas sites.bed --output results.csv

        scqpas --config custom.yaml --bam sample.bam --gtf annotation.gtf --pas sites.bed --output results.csv
    """

    try:
        # Load configuration
        try:
            config_manager = ConfigManager(config_path=config)
        except Exception as e:
            raise click.ClickException(f"Configuration error: {e}")

        # Apply CLI overrides to config values (CLI args take precedence)
        if percentage_threshold is None:
            percentage_threshold = config_manager.get(
                "polya_detection", "percentage_threshold"
            )

        if length_threshold is None:
            length_threshold = config_manager.get(
                "polya_detection", "length_threshold"
            )

        if terminal_exon_extension is None:
            terminal_exon_extension = config_manager.get(
                "annotation", "terminal_exon_extension"
            )

        if stringency is None:
            stringency = config_manager.get(
                "pas_filtering", "stringency"
            )

        if log_level is None:
            log_level = config_manager.get("logging", "default_level", "INFO")

        if log_file is None:
            log_file = config_manager.get("logging", "default_file", "scqpas.log")

        if output is None:
            output = config_manager.get("output", "default_output_file", "distances.csv")

        # Configure logging based on user's options
        log_level_int = getattr(logging, log_level.upper())
        # If log_file is empty string, pass None to disable file logging
        log_file_arg = log_file if log_file else None
        configure_logging(log_level=log_level_int, log_file=log_file_arg)

        run_pipeline(
            bam_path=bam,
            gtf_path=gtf,
            pas_bed_path=pas,
            output_path=output,
            percentage_threshold=percentage_threshold,
            length_threshold=length_threshold,
            terminal_exon_extension=terminal_exon_extension,
            stringency=stringency,
            region=region,
            config_manager=config_manager,
            debug_output_dir=debug_output,
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
