import logging
from typing import Optional, Union

import pandas as pd
import numpy as np

from .extract_cigar_BAM import extract_cigar_n_metrics

logger = logging.getLogger(__name__)


def _df_distances(
    reads_path: str,
    valid_reads_path: str,
    intronic_reads_path: str,
    out_distances_path: str,
) -> None:
    """
    Load read data and calculate distances from reads to cleavage sites.

    Parameters
    ----------
    reads_path : str
        Path to reads CSV file
    valid_reads_path : str
        Path to valid reads (filtered) BED format file
    intronic_reads_path : str
        Path to intronic reads intersection BED format file
    out_distances_path : str
        Path to output distances CSV file

    Returns
    -------
    None
        Writes results to out_distances_path
    """
    df = pd.read_csv(reads_path, index_col=0)
    # print(df)
    valid_reads = pd.read_csv(valid_reads_path, sep="\t")
    # print(valid_reads)
    intronic_reads = pd.read_csv(
        intronic_reads_path,
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
    # print(intronic_reads)
    # Determine the best CPA site
    df_cpa_sites = _get_cpa_sites(df)
    # print(df_cpa_sites)
    # best_cpa_site = best_cpa(df_cpa_sites)
    # print('Best CPA site: ' + str(best_cpa_site))
    df_distance = calculate_distances(valid_reads, intronic_reads)
    df_distance.to_csv(out_distances_path, index=False)


####### Modify so that it takes multiple CPA sites #######


def _get_cpa_sites(reads_df: pd.DataFrame) -> pd.DataFrame:
    """
    Count occurrences of each cleavage/polyadenylation site by chromosome and strand.

    Groups reads by chromosome, CPA site position, and strand to identify all putative
    cleavage sites and their supporting read counts. This allows distinguishing CPAs
    at the same genomic position on opposite strands (e.g., opposite gene directions).

    Parameters
    ----------
    reads_df : pd.DataFrame
        DataFrame containing read information with columns 'chr', 'cpa_site', and 'strand'.
        Typically the output from extract_reads() in extract_reads_BAM.py.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['chr', 'cpa_site', 'strand', 'supporting_read_count']
        showing the number of reads supporting each (chr, cpa_site, strand) tuple.
        Sorted by supporting_read_count in descending order.
    """
    # Remove rows with missing CPA sites
    df_valid = reads_df.dropna(subset=["cpa_site"])

    if df_valid.empty:
        raise RuntimeError("No valid CPA sites found in reads")

    # Count occurrences of each (chr, cpa_site, strand) tuple
    cpa_counts = (
        df_valid.groupby(["chr", "cpa_site", "strand"])
        .size()
        .reset_index(name="supporting_read_count")
    )

    # Sort by supporting read count (descending) for easy selection of best CPA
    cpa_counts = cpa_counts.sort_values(
        "supporting_read_count", ascending=False
    ).reset_index(drop=True)

    return cpa_counts


def _best_cpa(
    cpa_sites_df: pd.DataFrame,
) -> Optional[tuple]:
    """
    Identify the cleavage/polyadenylation site with the most supporting reads.

    Selects the (chr, cpa_site, strand) tuple with the highest supporting read count.
    This represents the most probable cleavage site in the dataset, including strand
    information to distinguish CPAs at the same position on opposite strands.

    Parameters
    ----------
    cpa_sites_df : pd.DataFrame
        DataFrame with columns ['chr', 'cpa_site', 'strand', 'supporting_read_count']
        as returned by get_cpa_sites(). Should already be sorted by
        supporting_read_count in descending order.

    Returns
    -------
    tuple or None
        (chr, cpa_site_position, strand) for the CPA site with highest read count,
        or None if no CPA sites are found
    """
    if cpa_sites_df.empty:
        return None

    # Get the first row (highest supporting read count due to sorting)
    best_row = cpa_sites_df.iloc[0]
    best_chr = best_row["chr"]
    best_pos = best_row["cpa_site"]
    best_strand = best_row["strand"]

    return (best_chr, best_pos, best_strand)


def calculate_distances(
    valid_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calculate distances from read start to polyadenylation site (PAS).

    Computes the distance from the 5' start of each read to the PAS position.
    Note: The input valid_reads should have PAS-adjusted coordinates where:
    - Forward strand (+): end_read = PAS position
    - Reverse strand (-): start_read = PAS position

    These adjusted coordinates are set by adjust_read_ends() and passed through
    the bedtools intersection results.

    Parameters
    ----------
    valid_reads : pd.DataFrame
        DataFrame with validated reads including columns:
        'read_id', 'transcript_id', 'pas_id', 'rs_id', 'intronic', 'CIGAR_N', 'strand'
        Note: start_read, end_read should be PAS-adjusted coordinates (from rs_pas_adj)
    intronic_reads : pd.DataFrame
        Intersection results with columns:
        'chr_read', 'start_read', 'end_read', 'read_id', 'pas_id', 'transcript_id', 'length_intron', etc.
    output_csv : str, optional
        Path to output CSV file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        DataFrame with distance calculations and supporting information.
        'distance' column = genomic distance from read start to PAS position,
        with intron lengths subtracted for intronic reads.
    """

    # Create a new DataFrame for storing distances
    distance_df = valid_reads[
        ["read_id", "rs_id", "transcript_id", "pas_id", "intronic", "CIGAR_N"]
    ].copy()
    # print(distance_df)

    # Merging to include total intron lengths in the distance_df
    distance_df = distance_df.merge(
        intronic_reads.groupby(["read_id", "transcript_id", "pas_id"])["length_intron"]
        .sum()
        .reset_index(),
        on=["read_id", "transcript_id", "pas_id"],
        how="left",
    )
    distance_df.rename(columns={"length_intron": "total_intron_len"}, inplace=True)

    # Calculate distance from read start to PAS position
    # After adjust_read_ends(), coordinates are set such that:
    # - For + strand: end_read = PAS position, start_read = original 5' start
    # - For - strand: start_read = PAS position - 1, end_read = original 5' extent (high coord)
    # Therefore, same formula works for both strands:
    distance_df["distance"] = (
        valid_reads["end_read"].values - valid_reads["start_read"].values
    )

    # For intronic reads, subtract total intron length from distance
    # (introns don't contribute to read-to-PAS distance)
    intronic_mask = valid_reads["intronic"] == True
    distance_df.loc[intronic_mask, "distance"] -= distance_df.loc[
        intronic_mask, "total_intron_len"
    ]

    # cols: "read_id", "rs_id", "transcript_id", "pas_id", "intronic", "CIGAR_N", "total_intron_len", "distance"

    # print(distance_df)

    if output_csv is not None:
        distance_df.to_csv(output_csv, index=False)

    return distance_df


def calculate_polyA_distances(
    polyA_reads_indexed: pd.DataFrame,
    cigar_metrics: pd.DataFrame,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calculate distances from polyA reads to their CPA sites (from atlas PAS).

    **Input structure:**
    - polyA_reads_indexed: one row per (read, cpa_site) combination
      Columns: read_id, rs_id, chr, start, end, strand, cpa_site, pas_id, UMI, CB

    This allows a single read to potentially have multiple cpa_sites if the RS contains
    reads with different cleavage sites.

    **Distance calculation:**
    - Forward strand (+): distance = cpa_site - read_start
      (distance from 5' read start to detected CPA position)
    - Reverse strand (-): distance = read_end - cpa_site
      (distance from 3' read end to detected CPA position, accounting for reverse direction)

    **CIGAR correction:**
    After initial distance calculation, subtract total CIGAR N length (gaps that don't
    contribute to read-to-CPA distance).

    **Output format:**
    Matches non-polyA valid_reads format for integration:
    - read_id, rs_id, transcript_id (NaN), pas_id, intronic (False), CIGAR_N,
      total_intron_len (= total_cigar_N_length), distance, is_polyA_RS (True)

    Parameters
    ----------
    polyA_reads_indexed : pd.DataFrame
        One row per (read, cpa_site) pair from filter_polyA_by_pas().
        Columns must include: read_id, rs_id, chr, start, end, strand, cpa_site, pas_id
    cigar_metrics : pd.DataFrame
        CIGAR metrics from extract_cigar_n_metrics() (in extract_cigar_BAM.py).
        Columns: read_id, CIGAR_N, total_cigar_N_length
    output_csv : str, optional
        Path to save output CSV. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        Distance calculations for polyA reads with columns:
        read_id, rs_id, transcript_id, pas_id, intronic, CIGAR_N, total_intron_len,
        distance, is_polyA_RS
    """

    # Start with the indexed reads
    result_df = polyA_reads_indexed[
        ["read_id", "rs_id", "start", "end", "strand", "cpa_site", "pas_id"]
    ].copy()

    # Merge CIGAR metrics
    result_df = result_df.merge(cigar_metrics, on="read_id", how="left")

    # Calculate distance based on strand
    # For + strand: cpa_site (3' end) - read_start (5' start)
    # For - strand: read_end (5' start on reverse) - cpa_site (3' end on reverse)
    mask_plus = result_df["strand"] == "+"
    mask_minus = result_df["strand"] == "-"

    result_df["distance"] = np.nan
    result_df.loc[mask_plus, "distance"] = (
        result_df.loc[mask_plus, "cpa_site"].values
        - result_df.loc[mask_plus, "start"].values
    )
    result_df.loc[mask_minus, "distance"] = (
        result_df.loc[mask_minus, "end"].values
        - result_df.loc[mask_minus, "cpa_site"].values
    )

    # Subtract CIGAR N lengths from distance
    # (gaps don't contribute to genomic distance from read to CPA)
    result_df["distance"] = result_df["distance"] - result_df["total_cigar_N_length"]

    # Prepare output with required columns
    result_df["transcript_id"] = np.nan  # polyA reads have no transcript assignment
    result_df["intronic"] = False  # polyA uses direct CPA, not transcript intron-based
    result_df["is_polyA_RS"] = True

    # Rename for consistency with other output
    result_df = result_df.rename(columns={"total_cigar_N_length": "total_intron_len"})

    # Select and order columns to match non-polyA format
    output_cols = [
        "read_id",
        "rs_id",
        "transcript_id",
        "pas_id",
        "intronic",
        "CIGAR_N",
        "total_intron_len",
        "distance",
        "is_polyA_RS",
    ]
    result_df = result_df[output_cols]

    if output_csv is not None:
        result_df.to_csv(output_csv, index=False)

    return result_df
