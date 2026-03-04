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
    required=False,
    default=None,
    help="Path to polyadenylation sites BED file (DEPRECATED - PAS now auto-detected from reads). Kept for backward compatibility.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default=None,
    help="Output file path (CSV format). Uses config default if not provided",
)
@click.option(
    "--percentage-threshold",
    type=int,
    default=None,
    help="Min %% of A nucleotides in polyA region (0-100). Uses config default if not provided",
)
@click.option(
    "--length-threshold",
    type=int,
    default=None,
    help="Min length of soft-clipped region at 3' end (bp). Uses config default if not provided",
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
def main(
    config: Optional[str],
    bam: str,
    gtf: str,
    pas: Optional[str],
    output: Optional[str],
    percentage_threshold: Optional[int],
    length_threshold: Optional[int],
    log_level: Optional[str],
    log_file: Optional[str],
) -> None:
    """
    Calculate distances from cleavage sites to polyadenylation sites.

    Processes a BAM file against a GTF annotation to quantify distances between
    reads and the automatically detected polyadenylation site.

    The PAS site is now automatically detected from the read data by identifying
    the cleavage position with the most supporting evidence. This data-driven
    approach eliminates the need for manual PAS specification.

    Note: The --pas argument is deprecated and no longer used. PAS is determined
    from read cleavage patterns instead.

    All processing is in-memory except temporary files for bedtools (auto-cleaned).

    Configuration can be provided via:
    - Custom YAML config file (--config)
    - Defaults from scqpas/config/defaults.yaml
    - Command-line arguments (override config)

    Example:

        scqpas --bam sample.bam --gtf annotation.gtf --output results.csv

        scqpas --config custom.yaml --bam sample.bam --gtf annotation.gtf --output results.csv
    """

    try:
        # Load configuration
        try:
            config_manager = ConfigManager(config_path=config)
        except Exception as e:
            raise click.ClickException(f"Configuration error: {e}")

        # Apply CLI overrides to config values (CLI args take precedence)
        if output is None:
            output = config_manager.get(
                "output", "default_output_file", "distances.csv"
            )

        if percentage_threshold is None:
            percentage_threshold = config_manager.get(
                "polya_detection", "percentage_threshold", 80
            )

        if length_threshold is None:
            length_threshold = config_manager.get(
                "polya_detection", "length_threshold", 5
            )

        if log_level is None:
            log_level = config_manager.get("logging", "default_level", "INFO")

        if log_file is None:
            log_file = config_manager.get("logging", "default_file", "scqpas.log")

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
            config_manager=config_manager,
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
