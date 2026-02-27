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

    unmatched_list = []

    # Handle both DataFrame and file path inputs
    if isinstance(df_input, str):
        df = pd.read_csv(df_input, sep=",")
    else:
        df = df_input.copy()

    for _, row in df.iterrows():
        cig_segments = re.findall(r"(\d+)([MIDNSHPX=])", row["CIGAR"])

        current_start = row["start"]
        for length, operation in cig_segments:
            length = int(length)
            if operation == "N":
                unmatched_row = {
                    "chr": row["chr"],
                    "start": current_start,
                    "end": current_start + length,
                    "read_id": row["read_id"],
                    "strand": row["strand"],
                    "CIGAR": row["CIGAR"],
                }
                unmatched_list.append(unmatched_row)
            elif operation in ["M", "D", "=", "X"]:
                current_start += length

    unmatched_df = pd.DataFrame(unmatched_list)

    if bed_output is not None:
        unmatched_df.to_csv(bed_output, sep="\t", header=False, index=False)

    return unmatched_df


def filter_by_cigar(
    reads: Union[pd.DataFrame, str],
    filtered_reads: Union[pd.DataFrame, str],
    intronic_reads: Union[pd.DataFrame, str],
    cigar_df: pd.DataFrame,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Filter reads by CIGAR string validation against annotation.

    - Keeps reads where CIGAR-derived introns match annotated introns
    - Keeps reads without introns (no N in CIGAR) within genes
    - Keeps reads close to cleavage site

    Parameters
    ----------
    reads : pd.DataFrame or str
        Reads DataFrame or path to reads CSV
    filtered_reads : pd.DataFrame or str
        Filtered reads (BED format) or path to BED file
    intronic_reads : pd.DataFrame or str
        Intersection of reads with introns (BED format) or path to BED file
    cigar_df : pd.DataFrame
        CIGAR-derived introns DataFrame
    output_csv : str, optional
        Path to output TSV file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        Valid reads with CIGAR validation columns added
    """

    # Handle both DataFrame and file path inputs
    if isinstance(reads, str):
        reads = pd.read_csv(reads, sep=",")
    else:
        reads = reads.copy()

    if isinstance(filtered_reads, str):
        filtered_reads = pd.read_csv(
            filtered_reads,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "read_id", "dummy", "strand"],
        )
    else:
        filtered_reads = filtered_reads.copy()
        # If DataFrame from bedtools with numeric column names (0,1,2,...), rename them
        if isinstance(filtered_reads.columns[0], int):
            filtered_reads.columns = [
                "chr",
                "start",
                "end",
                "read_id",
                "dummy",
                "strand",
            ]

    if isinstance(intronic_reads, str):
        intronic_reads = pd.read_csv(
            intronic_reads,
            sep="\t",
            header=None,
            names=[
                "chr_read",
                "start_read",
                "end_read",
                "read_id",
                "score_read",
                "strand_read",
                "chr_intron",
                "start_intron",
                "end_intron",
                "intron_id",
                "length_intron",
                "strand_intron",
                "gene_id",
            ],
        )
    else:
        intronic_reads = intronic_reads.copy()
        # If DataFrame from bedtools with numeric column names (0,1,2,...), rename them
        if isinstance(intronic_reads.columns[0], int):
            intronic_reads.columns = [
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

    cigar_df = cigar_df.copy()
    cigar_df.columns = ["chr_cf", "start_cf", "end_cf", "read_id", "strand_cf", "CIGAR"]
    intronic_reads["transcript_id"] = intronic_reads["intron_id"].str.split("_").str[0]

    # Merge DataFrames on read_id
    merged_df = pd.merge(
        cigar_df, intronic_reads, on="read_id", suffixes=("_cf", "_intron"), how="left"
    )
    # Create a boolean series for matches
    merged_df["match"] = (merged_df["start_cf"] == merged_df["start_intron"]) & (
        merged_df["end_cf"] == merged_df["end_intron"]
    )

    # Group by read_id and transcript_id to check if all introns matched
    intronic_N = (
        merged_df.groupby(["read_id", "transcript_id"])
        .agg(all_match=("match", "all"))
        .reset_index()
        .query("all_match")
        .drop(columns="all_match")
    )

    # Keep reads from df that do not have an N in their CIGAR string (i.e., no introns)
    cigar_noN = reads[reads["CIGAR"].str.contains("N") == False]["read_id"]
    intronic_noN = (
        intronic_reads[intronic_reads["read_id"].isin(cigar_noN)][
            ["read_id", "transcript_id"]
        ]
        .groupby(["read_id", "transcript_id"], as_index=False)
        .first()
    )

    # Mark them as intronic
    intronic_N["intronic"] = True
    intronic_N["CIGAR_N"] = True
    intronic_noN["intronic"] = True
    intronic_noN["CIGAR_N"] = False

    # Handle reads close to the cleavage site -> excluding those already investigated as intronic
    excluded_reads = pd.concat([intronic_reads[["read_id"]], cigar_df[["read_id"]]])
    close_reads_list = filtered_reads["read_id"]
    close_reads = close_reads_list[~close_reads_list.isin(excluded_reads["read_id"])]

    close_reads = pd.DataFrame(
        {
            "read_id": close_reads,
            "transcript_id": np.nan,
            "intronic": False,
            "CIGAR_N": False,
        }
    )

    # Combine intronic reads with valid CIGAR N, reads without CIGAR N, and non-intronic reads
    final_reads = pd.concat([intronic_N, intronic_noN, close_reads], ignore_index=True)
    # Merge final_reads with filtered_reads to keep chr, start, end
    final_reads = pd.merge(
        final_reads,
        filtered_reads[["read_id", "chr", "start", "end", "strand"]],
        on="read_id",
        how="left",
    )

    if output_csv is not None:
        final_reads.to_csv(output_csv, sep="\t", index=False)

    return final_reads
