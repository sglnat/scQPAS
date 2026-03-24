"""
Core pipeline logic for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module contains the main pipeline orchestration logic, independent of CLI interface.
All processing is in-memory except temporary files required for bedtools (auto-cleaned).
"""

import gc
import logging
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

from .config_manager import ConfigManager

from .extract_reads_BAM import extract_reads
from .extract_annotation_GTF import (
    load_gtf,
    extract_genes,
    extract_transcripts,
    extract_exons,
    calculate_introns,
    reads_to_bed,
)
from .extract_cigar_BAM import get_cigar_bed, filter_by_cigar, extract_cigar_n_metrics
from .calculate_distances import calculate_distances, calculate_polyA_distances
from .bedtools_intersections import run_bedtools_intersect
from .process_pas_BED import (
    add_read_set_ids,
    filter_by_transcript_association,
    get_adj_gtf,
    load_pas,
    assign_rs_pas,
    adjust_read_ends,
    filter_polyA_by_pas,
)

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: str,
    gtf_path: str,
    pas_bed_path: str,
    output_path: Optional[str] = None,
    percentage_threshold: Optional[int] = None,
    length_threshold: Optional[int] = None,
    terminal_exon_extension: Optional[int] = None,
    stringency: Optional[int] = None,
    region: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
    debug_output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Execute the complete scQPAS pipeline.

    Processes a BAM file against a GTF annotation in 9 steps:
    1. Extract polyA-containing reads from BAM
    2. Extract annotation (genes, exons, introns) from GTF
    3. Convert reads to BED format
    4. Extract CIGAR-derived introns from reads
    5. Intersect reads with annotated transcripts
    6. Load and prepare PAS data (shared by both pathways)
    7. PATHWAY 1 - Non-polyA reads: assign PAS, intersect introns, calculate distances (if transcript intersections exist)
    8. PATHWAY 2 - PolyA reads: filter to PAS, extract CIGAR metrics, calculate distances (independent)
    9. Combine results from both pathways

    The pipeline processes non-polyA and polyA reads separately:
    - Non-polyA: uses transcript context, intron annotations, and detected CPA sites
    - PolyA: uses direct PAS matching from sc PolyASite Atlas and CIGAR-derived metrics

    Parameters
    ----------
    bam_path : str
        Path to input BAM file (must be indexed)
    gtf_path : str
        Path to GTF annotation file
    pas_bed_path : str
        Path to polyadenylation sites BED file
    output_path : str, optional
        Output file path for results CSV. If None, uses "distances.csv" from config
    percentage_threshold : int, optional
        Min % of A nucleotides in polyA region (0-100).
        If None, retrieved from config. Default config: 80
    length_threshold : int, optional
        Min length of soft-clipped region at 3' end (bp).
        If None, retrieved from config. Default config: 5
    terminal_exon_extension : int, optional
        Length in base pairs to extend terminal exons for PAS capture.
        If None, retrieved from config. Default config: 1000
    stringency : int, optional
        Minimum stringency level for PAS from sc PolyASite Atlas.
        Only PAS with stringency >= this value will be retained.
        If None, retrieved from config. Default config: 80
    region : str, optional
        Optional genomic region for PAS filtering in format "chr:start-end" (e.g., "chr12:15671-26783")
        or just chromosome (e.g., "chr12"). If None, no coordinate filtering applied.
        Default: None (loads entire PAS file)
    config_manager : ConfigManager, optional
        ConfigManager instance for accessing pipeline configuration. Default: None
    debug_output_dir : str, optional
        Directory to save intermediate pipeline files for debugging/inspection.
        If None, intermediate files are not saved (in-memory only).
        All intermediate DataFrames will be written as CSV files in this directory.

    Returns
    -------
    pd.DataFrame
        Distance DataFrame with calculated distances from read start to PAS

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

    # Get thresholds from config if not provided
    if percentage_threshold is None:
        if config_manager is None:
            raise ValueError(
                "ConfigManager is required when percentage_threshold is None. "
                "Either provide percentage_threshold explicitly or pass config_manager."
            )
        percentage_threshold = config_manager.get(
            "polya_detection", "percentage_threshold"
        )

    if length_threshold is None:
        if config_manager is None:
            raise ValueError(
                "ConfigManager is required when length_threshold is None. "
                "Either provide length_threshold explicitly or pass config_manager."
            )
        length_threshold = config_manager.get("polya_detection", "length_threshold")

    if stringency is None:
        if config_manager is None:
            raise ValueError(
                "ConfigManager is required when stringency is None. "
                "Either provide stringency explicitly or pass config_manager."
            )
        stringency = config_manager.get("pas_filtering", "stringency")

    if output_path is None:
        if config_manager is None:
            config_manager = ConfigManager()
        output_path = config_manager.get(
            "output", "default_output_file", "distances.csv"
        )

    # Create debug output directory if specified
    if debug_output_dir is not None:
        Path(debug_output_dir).mkdir(parents=True, exist_ok=True)
        logger.info(f"Debug output directory: {Path(debug_output_dir).absolute()}")

    with tempfile.TemporaryDirectory() as tmpdir:

        # ========== STEP 1: EXTRACT READS ==========
        logger.info("[1/9] Extracting reads from BAM file...")
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            all_reads_df = extract_reads(
                sam,
                percentage_threshold=percentage_threshold,
                length_threshold=length_threshold,
                config_manager=config_manager,
                reads_output=(
                    f"{debug_output_dir}/01_reads.csv" if debug_output_dir else None
                ),
            )

        # Extract polyA reads and non-polyA reads into separate dataframes
        polyA_reads_df = all_reads_df[all_reads_df["is_polyA_RS"]]
        reads_df = all_reads_df[~all_reads_df["is_polyA_RS"]]

        polyA_count = len(polyA_reads_df)
        logger.info(f"      ✓ Extracted {len(all_reads_df)} total reads")
        logger.info(f"      ✓ Found {polyA_count} reads in polyA-containing read sets")

        # ========== STEP 2: EXTRACT ANNOTATION ==========
        logger.info("[2/9] Extracting annotation from GTF...")

        # Load GTF once for all extraction functions
        gtf_df = load_gtf(gtf_path)

        # Extract annotation from loaded GTF DataFrame
        # genes_df = extract_genes(gtf_df, config_manager=config_manager)
        exons_df = extract_exons(
            gtf_df,
            exons_output=(
                f"{debug_output_dir}/02a_exons.bed" if debug_output_dir else None
            ),
        )
        logger.info(f"      ✓ Extracted {len(exons_df)} exons")
        introns_df = calculate_introns(
            exons_df,
            introns_output=(
                f"{debug_output_dir}/02b_introns.bed" if debug_output_dir else None
            ),
        )
        logger.info(f"      ✓ Calculated {len(introns_df)} introns")
        transcripts_df = extract_transcripts(
            gtf_df,
            transcripts_output=(
                f"{debug_output_dir}/02c_transcripts.bed" if debug_output_dir else None
            ),
        )
        logger.info(f"      ✓ Extracted {len(transcripts_df)} transcripts")

        # logger.info(f"      ✓ Extracted {len(genes_df)} genes")
        # logger.info(f"      ✓ Using detected PAS: {pas_chr}:{pas_pos} ({pas_strand})")

        # ========== STEP 3: CONVERT READS TO BED ==========

        logger.info("[3/9] Converting reads to BED...")

        # Convert reads to BED format for bedtools intersections
        reads_bed = reads_to_bed(reads_df)

        # ========== STEP 4: EXTRACT CIGAR-DERIVED INTRONS ==========
        logger.info("[4/9] Parsing CIGAR strings...")

        # Pass DataFrame directly
        cigar_df = get_cigar_bed(reads_df)
        logger.info(f"      ✓ Extracted {len(cigar_df)} CIGAR-derived introns")

        # ========== STEP 5: BEDTOOLS INTERSECTIONS: READS-TRANSCRIPTS ==========
        logger.info("[5/9] Intersecting reads with annotated transcripts...")

        # Intersect reads with transcripts
        reads_transcripts_df = run_bedtools_intersect(
            reads_bed,
            transcripts_df,
            tmpdir,
            name_a="reads",
            name_b="transcripts",
            flags=["-wa", "-wb", "-s", "-loj", "-f", "1.0"],
            # only keep reads that are 100% within a transcript
            # keep reads without overlap (loj) to identify all RS that don't overlap with transcripts
            output_bed=(
                f"{debug_output_dir}/05a_reads_transcripts_intersect.bed"
                if debug_output_dir
                else None
            ),
        )

        # Assign column names: -wa outputs 6 cols, -wb outputs 7 cols
        if not reads_transcripts_df.empty:
            # Keep only the reads cols, transcript_id and gene_id
            reads_transcripts_df = reads_transcripts_df[
                [0, 1, 2, 3, 4, 5, 9, 12]
            ]  # Column selection creates new DF
            reads_transcripts_df.columns = [
                "chr_read",
                "start_read",
                "end_read",
                "read_id",
                "dummy",
                "strand",
                "transcript_id",
                "gene_id",
            ]

            # Add read set IDs (rs_id = CB_UMI) by tracing back to original reads
            reads_transcripts_df = add_read_set_ids(reads_transcripts_df, reads_df)

            # Filter reads based on transcript association within read sets:
            # Keep only RS-transcript combinations where all reads from RS are within the transcript
            reads_transcripts_df = filter_by_transcript_association(
                reads_transcripts_df
            )

            if debug_output_dir is not None:
                reads_transcripts_df.to_csv(
                    f"{debug_output_dir}/05b_reads_transcripts_filtered.csv",
                    index=False,
                )

            logger.info(
                f"      ✓ Found {len(reads_transcripts_df)} read-transcript intersections"
            )

        # ========== STEP 6: LOAD AND PREPARE PAS (SHARED BY BOTH PATHWAYS) ==========
        logger.info("[6/9] Loading and preparing PAS data...")

        # Transcripts with terminal exons extended to capture more PAS
        gtf_adj = get_adj_gtf(
            gtf_path,
            output_file=(
                f"{debug_output_dir}/06_gtf_adjusted.gtf"
                if debug_output_dir
                else None
            ),
            extension_length=terminal_exon_extension,
            config_manager=config_manager,
        )
        transcripts_adj = extract_transcripts(
            gtf_adj,
            transcripts_output=(
                f"{debug_output_dir}/06_transcripts_adjusted.bed"
                if debug_output_dir
                else None
            ),
        )
        logger.info(
            f"      ✓ Extracted {len(transcripts_adj)} adjusted transcripts (terminal exons extended by {terminal_exon_extension}bp)"
        )

        # Load PAS bed file from sc PolyASite Atlas v3.0
        pas_df = load_pas(pas_bed_path, stringency, region=region)
        logger.info(f"      ✓ Loaded {len(pas_df)} PAS from atlas")
        if region is not None:
            logger.info(f"      ✓ Filtered to {len(pas_df)} PAS in region {region}")

        # Intersect PAS with adjusted transcripts to assign PAS to transcripts
        pas_transcript_df = run_bedtools_intersect(
            pas_df,
            transcripts_adj,
            tmpdir,
            name_a="pas",
            name_b="transcripts_adj",
            flags=["-wa", "-wb", "-s"],
            output_bed=(
                f"{debug_output_dir}/06a_pas_transcripts_intersect.bed"
                if debug_output_dir
                else None
            ),
        )
        logger.info(
            f"      ✓ Found {len(pas_transcript_df)} PAS-transcript intersections"
        )

        pas_transcript_df = pas_transcript_df[
            [0, 1, 2, 3, 4, 5, 9, 12]
        ]  # Column selection creates new DF
        pas_transcript_df.columns = [
            "chr_pas",
            "start_pas",
            "end_pas",
            "pas_id",
            "gex",
            "strand_pas",
            "transcript_id",
            "gene_id",
        ]

        # ========== PATHWAY 1: NON-POLYA READS (CONDITIONAL ON TRANSCRIPT INTERSECTIONS) ==========
        distances_df = pd.DataFrame()
        if not reads_transcripts_df.empty:
            logger.info("[7/9] Processing non-polyA reads...")

            # Assign PAS to reads based on their associated transcripts and read sets
            rs_pas_df = assign_rs_pas(
                reads_transcripts_df,
                pas_transcript_df,
                config_manager=config_manager,
                output_csv=(
                    f"{debug_output_dir}/07a_rs_pas_assignments.csv"
                    if debug_output_dir
                    else None
                ),
            )
            logger.info(f"      ✓ Assigned PAS to {len(rs_pas_df)} read-set records")

            rs_pas_adj = adjust_read_ends(
                rs_pas_df,
                output_csv=(
                    f"{debug_output_dir}/07b_rs_pas_adjusted.csv"
                    if debug_output_dir
                    else None
                ),
            )
            logger.info(
                f"      ✓ Adjusted read coordinates for {len(rs_pas_adj)} read-set records"
            )

            # Intersect reads with introns
            reads_introns_df = run_bedtools_intersect(
                rs_pas_adj,
                introns_df,
                tmpdir,
                name_a="reads_pas",
                name_b="introns",
                flags=["-s", "-f", "1.0", "-wa", "-wb"],
                output_bed=(
                    f"{debug_output_dir}/07c_reads_introns_intersect.bed"
                    if debug_output_dir
                    else None
                ),
            )

            # Assign column names: -wa outputs A columns, -wb outputs B columns
            if not reads_introns_df.empty:
                logger.info(
                    f"      ✓ Found {len(reads_introns_df)} raw read-intron intersections"
                )

                # rs_pas_adj columns: chr_read, start_read, end_read, read_id, dummy, strand, transcript_id, pas_id, rs_id, chr_pas, pos_pas, strand_pas
                # introns_df columns: chr, start, end, intron_id, length_intron, strand, gene_id
                reads_introns_df.columns = [
                    "chr_read",
                    "start_read",
                    "end_read",
                    "read_id",
                    "dummy",
                    "strand",
                    "transcript_id",
                    "pas_id",
                    "rs_id",
                    "chr_pas",
                    "pos_pas",
                    "strand_pas",
                    "chr_intron",
                    "start_intron",
                    "end_intron",
                    "intron_id",
                    "length_intron",
                    "strand_intron",
                    "gene_id",
                ]

                reads_introns_df["transcript_id_pas"] = (
                    reads_introns_df["intron_id"].str.split("_").str[0]
                )
                # reads_introns_df['transcript_id'] = reads_introns_df['transcript_id'].astype(str)
                reads_introns_df = reads_introns_df[
                    reads_introns_df["transcript_id"]
                    == reads_introns_df["transcript_id_pas"]
                ]
                logger.info(
                    f"      ✓ Filtered to {len(reads_introns_df)} matching transcript-intron pairs"
                )

                reads_introns_df.drop(columns=["transcript_id_pas"], inplace=True)

                # Drop unnecessary columns from bedtools output to save memory
                reads_introns_df.drop(
                    columns=["dummy", "chr_pas", "pos_pas", "strand_pas", "gene_id"],
                    inplace=True,
                )
                # Now has only: chr_read, start_read, end_read, read_id, strand, transcript_id, pas_id, rs_id, chr_intron, start_intron, end_intron, intron_id, length_intron, strand_intron

                if debug_output_dir is not None:
                    reads_introns_df.to_csv(
                        f"{debug_output_dir}/07d_reads_introns_filtered.csv",
                        index=False,
                    )

                logger.info(
                    f"      ✓ Retained final {len(reads_introns_df)} read-intron intersections"
                )

                # Free memory before filter_by_cigar
                del gtf_df, exons_df, introns_df, transcripts_df, reads_bed
                del transcripts_adj, pas_transcript_df
                del reads_transcripts_df, rs_pas_df
                gc.collect()

                valid_reads = filter_by_cigar(
                    reads_df, rs_pas_adj, reads_introns_df, cigar_df
                )
                logger.info(
                    f"      ✓ Validated {len(valid_reads)} reads against CIGAR introns"
                )
                # cols: "read_id", "transcript_id", "pas_id", "rs_id", "intronic", "CIGAR_N", "chr_read", "start_read", "end_read", "strand"

                if debug_output_dir is not None:
                    valid_reads.to_csv(
                        f"{debug_output_dir}/07e_valid_reads_filtered.csv", index=False
                    )

                n_intronic = (
                    valid_reads["intronic"].sum()
                    if "intronic" in valid_reads.columns
                    else 0
                )
                logger.info(
                    f"      ✓ Retained {len(valid_reads)} valid reads ({n_intronic} intronic, {len(valid_reads) - n_intronic} non-intronic)"
                )

                distances_df = calculate_distances(
                    valid_reads,
                    reads_introns_df,
                    output_csv=(
                        f"{debug_output_dir}/07f_distances_non_polyA.csv"
                        if debug_output_dir
                        else None
                    ),
                )
                logger.info(
                    f"      ✓ Calculated distances for {len(distances_df)} non-polyA reads"
                )

                # Free up memory
                del valid_reads
                gc.collect()

            else:
                # No intron intersections found, start with empty distances df
                logger.warning("      ✗ No read-intron intersections found")
                distances_df = pd.DataFrame()
        else:
            # No read-transcript intersections found, start with empty distances df
            logger.warning("      ✗ No read-transcript intersections found")
            distances_df = pd.DataFrame()

        # ========== PATHWAY 2: POLYA READS ==========
        logger.info("[8/9] Processing polyA reads and calculating distances...")

        polyA_distances_df = None
        if not polyA_reads_df.empty:
            # Filter polyA reads by atlas PAS
            polyA_reads_indexed = filter_polyA_by_pas(
                polyA_reads_df,
                pas_df,
                output_csv=(
                    f"{debug_output_dir}/08a_polyA_reads_indexed.csv"
                    if debug_output_dir
                    else None
                ),
            )
            logger.info(f"      ✓ Indexed {len(polyA_reads_indexed)} polyA reads to PAS")

            # Extract CIGAR metrics for polyA reads
            polyA_cigar_metrics = extract_cigar_n_metrics(polyA_reads_df)

            # Calculate distances from polyA reads to CPA sites
            polyA_distances_df = calculate_polyA_distances(
                polyA_reads_indexed,
                polyA_cigar_metrics,
                output_csv=(
                    f"{debug_output_dir}/08b_polyA_distances.csv"
                    if debug_output_dir
                    else None
                ),
            )
            logger.info(
                f"      ✓ Calculated distances for {len(polyA_distances_df)} polyA reads"
            )
            del polyA_reads_indexed, polyA_cigar_metrics
            gc.collect()
        else:
            logger.info("      ✓ No polyA reads to process")

        # ========== FINAL MERGE: COMBINE NON-POLYA AND POLYA RESULTS ==========
        logger.info("[9/9] Combining results...")

        # Merge based on what distance dataframes we have
        if not distances_df.empty and polyA_distances_df is not None and not polyA_distances_df.empty:
            # Both pathways produced results
            distances_df = pd.concat([distances_df, polyA_distances_df], ignore_index=True)
            logger.info(
                f"      ✓ Combined {len(distances_df)} total distances (non-polyA + polyA)"
            )
        elif polyA_distances_df is not None and not polyA_distances_df.empty:
            # Only polyA pathway produced results
            distances_df = polyA_distances_df
            logger.info(f"      ✓ Using {len(distances_df)} polyA distances only")
        # else: distances_df already has non-polyA results or is empty (handled above)

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
