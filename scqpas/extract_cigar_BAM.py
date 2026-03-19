import logging
from typing import Optional, Union

import pandas as pd
import numpy as np
import re

logger = logging.getLogger(__name__)


def get_cigar_bed(df_input: Union[pd.DataFrame, str], bed_output: Optional[str] = None) -> pd.DataFrame:
    """
    Extract CIGAR-derived introns (N operations) from reads.

    Parameters
    ----------
    df_input : pd.DataFrame or str
        Reads DataFrame or path to reads CSV file.
        Must have columns: 'chr', 'start', 'read_id', 'strand', 'CIGAR'
    bed_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        CIGAR-derived introns with columns: chr, start, end, read_id, strand, CIGAR
    """

    # Handle both DataFrame and file path inputs
    if isinstance(df_input, str):
        df = pd.read_csv(df_input, sep=",")
    else:
        df = df_input

    # Vectorized CIGAR parsing using apply instead of iterrows
    def parse_cigar_row(row):
        """Extract N introns from a single CIGAR string"""
        introns = []
        cig_segments = re.findall(r"(\d+)([MIDNSHPX=])", row["CIGAR"])
        current_start = row["start"]
        
        for length, operation in cig_segments:
            length = int(length)
            if operation == "N":
                introns.append({
                    "chr": row["chr"],
                    "start": current_start,
                    "end": current_start + length,
                    "read_id": row["read_id"],
                    "strand": row["strand"],
                    "CIGAR": row["CIGAR"],
                })
            elif operation in ["M", "D", "=", "X"]:
                current_start += length
        
        return introns

    # Apply vectorized parsing - much faster than iterrows
    all_introns = []
    for introns in df.apply(parse_cigar_row, axis=1):
        all_introns.extend(introns)
    
    unmatched_df = pd.DataFrame(all_introns) if all_introns else pd.DataFrame()

    if bed_output is not None:
        unmatched_df.to_csv(bed_output, sep="\t", header=False, index=False)

    return unmatched_df


def filter_by_cigar(
    reads: pd.DataFrame,
    filtered_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame,
    cigar_df: pd.DataFrame,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Filter reads by CIGAR string validation against annotation.

    Validates that CIGAR-derived introns match annotated introns for each
    (read_id, transcript_id, rs_id) combination. Rejects entire (rs_id, transcript_id)
    pairs if ANY read in that set has mismatched introns.

    Parameters
    ----------
    reads : pd.DataFrame
        Reads DataFrame with columns: read_id, chr, start, end, strand, CIGAR
    filtered_reads : pd.DataFrame
        Read-transcript-PAS assignments with columns: 
        chr_read, start_read, end_read, read_id, strand, transcript_id, pas_id, rs_id
    intronic_reads : pd.DataFrame
        Intersection of reads with annotated introns with columns:
        chr_read, start_read, end_read, read_id, strand, transcript_id, pas_id, rs_id,
        chr_intron, start_intron, end_intron, intron_id, length_intron, strand_intron
    cigar_df : pd.DataFrame
        CIGAR-derived introns (from get_cigar_bed) with columns:
        chr, start, end, read_id, strand, CIGAR
    output_csv : str, optional
        Path to output TSV file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        Valid reads with columns: read_id, transcript_id, pas_id, rs_id,
        chr_read, start_read, end_read, strand, intronic, CIGAR_N
    """

    # ========== STEP 1: SPLIT READS BY CIGAR N CONTENT ==========
    # Separate reads into two groups:
    # - reads_with_N: Have splice junctions (need CIGAR validation)
    # - reads_without_N: Continuous reads (no validation needed)
    
    has_N_mask = reads["CIGAR"].str.contains("N", na=False)
    reads_with_N_set = set(reads.loc[has_N_mask, "read_id"])
    reads_without_N_set = set(reads.loc[~has_N_mask, "read_id"])
    
    # Split intronic_reads by CIGAR N content
    intronic_noN = intronic_reads[intronic_reads["read_id"].isin(reads_without_N_set)].copy()
    intronic_reads_with_N = intronic_reads[intronic_reads["read_id"].isin(reads_with_N_set)].copy()

    # ========== STEP 2: BUILD CIGAR INTRON DICTIONARY (by read_id) ==========
    # Extract CIGAR-derived introns grouped by read_id
    # For each read, store a SET of (start, end) coordinate tuples
    # This allows O(1) lookup of whether a coordinate pair exists
    
    cigar_by_read = {}
    for read_id, group in cigar_df.groupby("read_id"):
        cigar_by_read[read_id] = set(zip(group["start"], group["end"]))
    
    # ========== STEP 3: BUILD ANNOTATED INTRON DICTIONARY (by read_id, transcript_id) ==========
    # For CIGAR validation, we need to check: does this read's CIGAR match the annotation
    # for this specific transcript?
    # DEDUPLICATION: Group by (read_id, transcript_id) only, not rs_id
    # Same read+transcript combination has identical CIGAR introns regardless of PAS assignment
    # Multiple PAS for same read+transcript will merge into one validation check
    # Store as: {(read_id, transcript_id) → set((start, end))}
    
    annotated_by_combo = {}
    for (read_id, transcript_id), group in intronic_reads_with_N.groupby(["read_id", "transcript_id"]):
        # Extract unique intron coordinates from annotation for this combo
        # Multiple pas_ids will merge into single intron set (they're all the same anyway)
        # Store as set of (start, end) tuples
        annotated_by_combo[(read_id, transcript_id)] = set(zip(group["start_intron"], group["end_intron"]))
    
    # Create mapping from (read_id, transcript_id) to rs_id for rejection tracking
    # Every (read_id, transcript_id) pair maps to exactly one rs_id (read sets don't split)
    # Use .first() to get first occurrence (all are identical for that pair)
    read_tx_to_rs = intronic_reads_with_N.groupby(["read_id", "transcript_id"])["rs_id"].first().to_dict()

    # ========== STEP 4: VALIDATE CIGAR INTRONS MATCH ANNOTATION ==========
    # For each (read_id, transcript_id) combo (deduplicated):
    # Check if ALL CIGAR introns from this read are present in the annotated introns
    # Validation rule: cigar_introns ⊆ annotated_introns (subset check)
    # If validation FAILS → mark entire (rs_id, transcript_id) as invalid
    # This ensures all reads with same UMI+CB (rs_id) are consistent for each transcript
    
    rejected_rs_transcript = set()  # Track (rs_id, transcript_id) pairs that fail validation
    
    for (read_id, transcript_id), annotated in annotated_by_combo.items():
        cigar_introns = cigar_by_read.get(read_id, set())
        # Check: do all CIGAR introns exist in the annotation for this (read, transcript)?
        if not cigar_introns.issubset(annotated):
            # VALIDATION FAILED: This read's CIGAR doesn't match the annotation
            # Look up rs_id for this (read_id, transcript_id) pair and reject both
            rs_id = read_tx_to_rs[(read_id, transcript_id)]
            # Remove the entire (rs_id, transcript_id) pair from consideration
            rejected_rs_transcript.add((rs_id, transcript_id))
    
    # ========== STEP 5: FILTER OUT INVALID (rs_id, transcript_id) PAIRS ==========
    # Remove all rows where (rs_id, transcript_id) pair was rejected
    # This means ALL reads from this read set (rs_id) are removed for this transcript
    
    if rejected_rs_transcript:
        # Convert rejected set to DataFrame for efficient merge
        rejected_df = pd.DataFrame(
            list(rejected_rs_transcript),
            columns=["rs_id", "transcript_id"]
        )
        # Left merge with indicator: identifies which rows are in rejected_df
        merged = intronic_reads_with_N.merge(
            rejected_df,
            on=["rs_id", "transcript_id"],
            how="left",
            indicator=True
        )
        # Keep only rows NOT in rejected_df (left_only = in intronic_reads_with_N but not in rejected_df)
        intronic_reads_with_N = intronic_reads_with_N[merged["_merge"] == "left_only"].copy()
        del merged  # Free memory

    # ========== STEP 6: PREPARE intronic_N (READS WITH N THAT PASSED VALIDATION) ==========
    # Extract unique (read_id, transcript_id, pas_id, rs_id) from validated reads
    # Add flags: intronic=True (these reads overlap introns), CIGAR_N=True (they have splicing)
    
    if not intronic_reads_with_N.empty:
        intronic_N = intronic_reads_with_N[["read_id", "transcript_id", "pas_id", "rs_id"]].drop_duplicates()
        intronic_N["intronic"] = True
        intronic_N["CIGAR_N"] = True
    else:
        intronic_N = pd.DataFrame(columns=["read_id", "transcript_id", "pas_id", "rs_id", "intronic", "CIGAR_N"])

    # ========== STEP 7: PREPARE intronic_noN (READS WITHOUT N) ==========
    # These reads have no splice junctions, so no CIGAR validation needed
    # They automatically pass because they're in intronic_reads (overlapping annotation)
    # Flag: intronic=True (overlap with introns), CIGAR_N=False (no splicing)
    
    if not intronic_noN.empty:
        intronic_noN = intronic_noN[["read_id", "transcript_id", "pas_id", "rs_id"]].drop_duplicates()
        intronic_noN["intronic"] = True
        intronic_noN["CIGAR_N"] = False
    else:
        intronic_noN = pd.DataFrame(columns=["read_id", "transcript_id", "pas_id", "rs_id", "intronic", "CIGAR_N"])

    # ========== STEP 8: HANDLE READS CLOSE TO CLEAVAGE SITE ==========
    # Find reads in filtered_reads that:
    # 1. Are NOT already in intronic_reads (not counted as intronic)
    # 2. Have (rs_id, transcript_id) that was NOT rejected in Step 5
    # These reads are close to the PAS but don't overlap introns
    
    close_reads_df = pd.DataFrame()
    
    if filtered_reads is not None and not filtered_reads.empty:
        # Build DataFrame of (read_id, rs_id, transcript_id) from intronic reads for efficient anti-join
        intronic_combos_df = pd.DataFrame()
        if not intronic_reads_with_N.empty:
            intronic_combos_df = pd.concat([
                intronic_combos_df,
                intronic_reads_with_N[["read_id", "rs_id", "transcript_id"]]
            ], ignore_index=True)
        if not intronic_noN.empty:
            intronic_combos_df = pd.concat([
                intronic_combos_df,
                intronic_noN[["read_id", "rs_id", "transcript_id"]]
            ], ignore_index=True)
        
        if not intronic_combos_df.empty:
            intronic_combos_df["_intronic"] = True
        
        # Left merge: identifies reads that are in filtered_reads but NOT in intronic_reads
        merged = filtered_reads.merge(
            intronic_combos_df,
            on=["read_id", "rs_id", "transcript_id"],
            how="left"
        )
        # Keep rows where _intronic is NaN (not in intronic reads)
        not_intronic = merged[merged["_intronic"].isna()].index
        
        # Also exclude rows where (rs_id, transcript_id) was rejected
        if rejected_rs_transcript:
            rejected_df = pd.DataFrame(
                list(rejected_rs_transcript),
                columns=["rs_id", "transcript_id"]
            )
            merged2 = filtered_reads.loc[not_intronic].merge(
                rejected_df,
                on=["rs_id", "transcript_id"],
                how="left",
                indicator=True
            )
            not_rejected_mask = merged2["_merge"] == "left_only"
            close_reads_slice = filtered_reads.loc[not_intronic][not_rejected_mask]
        else:
            close_reads_slice = filtered_reads.loc[not_intronic]
        
        if not close_reads_slice.empty:
            close_reads_df = close_reads_slice[["read_id", "transcript_id", "pas_id", "rs_id"]].drop_duplicates()
            close_reads_df = close_reads_df.copy()  # Avoid SettingWithCopyWarning
            close_reads_df["intronic"] = False
            close_reads_df["CIGAR_N"] = False

    # ========== STEP 9: COMBINE AND RETURN ==========
    # Concatenate all valid reads: intronic_N, intronic_noN, close_reads
    # Merge with filtered_reads to get genomic coordinates
    
    final_reads = pd.concat([intronic_N, intronic_noN, close_reads_df], ignore_index=True)
    
    # Clean up large intermediate DataFrames no longer needed
    del cigar_by_read, annotated_by_combo, intronic_reads_with_N, intronic_noN, intronic_N, close_reads_df
    
    if not final_reads.empty:
        # Merge with filtered_reads to get chr_read, start_read, end_read, strand info
        final_reads = final_reads.merge(
            filtered_reads[["read_id", "chr_read", "start_read", "end_read", "strand"]],
            on="read_id",
            how="left",
        )
    
    if output_csv is not None:
        final_reads.to_csv(output_csv, sep="\t", index=False)

    return final_reads
