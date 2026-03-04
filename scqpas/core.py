"""
Core pipeline logic for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module contains the main pipeline orchestration logic, independent of CLI interface.
All processing is in-memory except temporary files required for bedtools (auto-cleaned).
"""

import logging
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

from .config_manager import ConfigManager

from .extract_reads_BAM import extract_reads
from .extract_annotation_GTF import (
    extract_genes,
    extract_exons,
    calculate_introns,
    get_bed_from_df,
)
from .extract_cigar_BAM import get_cigar_bed, filter_by_cigar
from .calculate_distances import calculate_distances, detect_best_cpa_from_reads
from .bedtools_intersections import run_bedtools_intersect

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: str,
    gtf_path: str,
    pas_bed_path: Optional[str] = None,
    output_path: str = None,
    percentage_threshold: int = 80,
    length_threshold: int = 5,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Execute the complete scQPAS pipeline.

    Processes a BAM file against a GTF annotation in 6 steps:
    1. Extract polyA-containing reads from BAM
    2. Detect best CPA site from read data
    3. Extract annotation (genes, exons, introns) from GTF
    4. Extract CIGAR-derived introns from reads
    5. Intersect reads with annotated features (bedtools)
    6. Calculate distances from cleavage to detected PAS

    The polyadenylation site is automatically detected from the reads themselves
    by identifying the cleavage site with the most supporting reads. This data-driven
    approach replaces manual specification of the PAS position.

    Parameters
    ----------
    bam_path : str
        Path to input BAM file (must be indexed)
    gtf_path : str
        Path to GTF annotation file
    pas_bed_path : str, optional
        Path to polyadenylation sites BED file (DEPRECATED).
        Currently unused; PAS is determined from read data instead.
        Kept for backward compatibility; will be used in future for atlas filtering.
        Default: None
    output_path : str
        Output file path for results CSV
    percentage_threshold : int, optional
        Min % of A nucleotides in polyA region (0-100). Default: 80
    length_threshold : int, optional
        Min length of soft-clipped region at 3' end (bp). Default: 5
    config_manager : ConfigManager, optional
        ConfigManager instance for accessing pipeline configuration. Default: None

    Returns
    -------
    pd.DataFrame
        Distance DataFrame with calculated distances from detected PAS

    Raises
    ------
    FileNotFoundError
        If input files don't exist
    RuntimeError
        If pipeline processing fails or no polyA reads found to detect CPA site
    """

    logger.info("=" * 70)
    logger.info("scQPAS: Single-cell Quantification of PolyAdenylation Sites")
    logger.info("=" * 70)

    if pas_bed_path is not None:
        logger.warning(
            "--pas argument is deprecated and no longer used. "
            "PAS is automatically detected from read data instead."
        )

    with tempfile.TemporaryDirectory() as tmpdir:

        # ========== STEP 1: EXTRACT READS ==========
        logger.info("[1/6] Extracting reads from BAM file...")
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            reads_df = extract_reads(
                sam,
                percentage_threshold=percentage_threshold,
                length_threshold=length_threshold,
                config_manager=config_manager,
            )

        polyA_count = reads_df["is_polyA"].sum()
        logger.info(f"      ✓ Extracted {len(reads_df)} total reads")
        logger.info(f"      ✓ Found {polyA_count} polyA-containing reads")

        # ========== STEP 2: DETECT BEST CPA SITE FROM READS ==========
        logger.info("[2/6] Detecting best CPA site from read data...")
        pas_chr, pas_pos, pas_strand = detect_best_cpa_from_reads(reads_df)

        # ========== STEP 3: EXTRACT ANNOTATION ==========
        logger.info("[3/6] Extracting annotation from GTF...")

        # All functions now support in-memory operation (output=None)
        genes_df = extract_genes(gtf_path, config_manager=config_manager)
        exons_df = extract_exons(gtf_path, config_manager=config_manager)
        introns_df = calculate_introns(exons_df, config_manager=config_manager)

        logger.info(f"      ✓ Extracted {len(genes_df)} genes")
        logger.info(f"      ✓ Extracted {len(exons_df)} exons")
        logger.info(f"      ✓ Calculated {len(introns_df)} introns")
        logger.info(f"      ✓ Using detected PAS: {pas_chr}:{pas_pos} ({pas_strand})")

        # ========== STEP 4: EXTRACT CIGAR-DERIVED INTRONS ==========
        logger.info("[4/6] Parsing CIGAR strings...")

        # Pass DataFrame directly (no file conversion needed)
        cigar_df = get_cigar_bed(reads_df)
        logger.info(f"      ✓ Extracted {len(cigar_df)} CIGAR-derived introns")

        # ========== STEP 5: RUN BEDTOOLS INTERSECTIONS ==========
        logger.info("[5/6] Intersecting reads with annotated features (using detected CPA)...")

        # Convert reads and annotation to BED format for this PAS
        reads_bed_df = get_bed_from_df(
            reads_df, pas_pos, chr=pas_chr, config_manager=config_manager
        )
        genes_bed_df = genes_df[["chr", "start", "end", "gene_id", "dummy", "strand"]]
        introns_bed_df = introns_df[
            ["chr", "start", "end", "intron_id", "length_intron", "strand"]
        ]

        # Intersect reads with genes
        reads_genes_df = run_bedtools_intersect(
            reads_bed_df,
            genes_bed_df,
            tmpdir,
            name_a="reads",
            name_b="genes",
            flags=["-wa", "-s", "-f", "1.0"],
        )

        # Assign column names: -wa outputs columns from reads_bed (6 cols)
        if not reads_genes_df.empty:
            reads_genes_df.columns = [
                "chr",
                "start",
                "end",
                "read_id",
                "dummy",
                "strand",
            ]

            logger.info(f"      ✓ Found {len(reads_genes_df)} read-gene intersections")

            # Intersect reads with introns
            reads_introns_df = run_bedtools_intersect(
                reads_genes_df,
                introns_bed_df,
                tmpdir,
                name_a="reads_genes",
                name_b="introns",
                flags=["-s", "-F", "1.0", "-wa", "-wb"],
            )

            # Assign column names: -wa outputs A columns, -wb outputs B columns
            if not reads_introns_df.empty:
                reads_introns_df.columns = [
                    "chr_read",
                    "start_read",
                    "end_read",
                    "read_id",
                    "dummy_read",
                    "strand_read",
                    "chr_intron",
                    "start_intron",
                    "end_intron",
                    "intron_id",
                    "length_intron",
                    "strand_intron",
                ]

                logger.info(f"      ✓ Found {len(reads_introns_df)} read-intron intersections")

                # ========== STEP 6: CALCULATE DISTANCES ==========
                logger.info("[6/6] Calculating distances from CPA sites...")

                valid_reads = filter_by_cigar(
                    reads_df, reads_genes_df, reads_introns_df, cigar_df
                )

                logger.info(f"      ✓ Retained {len(valid_reads)} valid reads")
                logger.info(f"      ✓ Intronic reads: {valid_reads['intronic'].sum()}")

                distances_df = calculate_distances(valid_reads, reads_introns_df)

                logger.info(f"      ✓ Calculated distances for {len(distances_df)} reads")

                # Add PAS coordinates (chr and position)
                distances_df["pas_chr"] = pas_chr
                distances_df["pas_pos"] = pas_pos
            else:
                distances_df = pd.DataFrame()
                logger.warning("      ✗ No read-intron intersections found")
        else:
            distances_df = pd.DataFrame()
            logger.warning("      ✗ No read-gene intersections found")

    # ========== WRITE OUTPUT ==========
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    distances_df.to_csv(output_path_obj, index=False)

    logger.info("=" * 70)
    logger.info("✓ PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 70)
    logger.info(f"Results written to: {output_path_obj.absolute()}")
    logger.info(f"Total reads processed: {len(distances_df)}")
    logger.info("=" * 70)

    return distances_df

