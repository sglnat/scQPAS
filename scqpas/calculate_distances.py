import logging
from typing import Optional, Union

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def df_distances(
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
    df_cpa_sites = get_cpa_sites(df)
    # print(df_cpa_sites)
    # best_cpa_site = best_cpa(df_cpa_sites)
    # print('Best CPA site: ' + str(best_cpa_site))
    df_distance = calculate_distances(valid_reads, intronic_reads)
    df_distance.to_csv(out_distances_path, index=False)


####### Modify so that it takes multiple CPA sites #######


def get_cpa_sites(reads_df: pd.DataFrame) -> pd.DataFrame:
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


def best_cpa(
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
    valid_reads: pd.DataFrame, intronic_reads: pd.DataFrame, output_csv: Optional[str] = None
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
    distance_df = valid_reads[["read_id", "rs_id", "transcript_id", "pas_id", "intronic", "CIGAR_N"]].copy()
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

    # Calculate distance based on strand information
    distance_df["distance"] = None

    # Forward strand (+): distance = PAS position (end_read) - read start (start_read)
    mask_plus = valid_reads["strand"] == "+"
    distance_df.loc[mask_plus, "distance"] = (
        valid_reads.loc[mask_plus, "end_read"].values
        - valid_reads.loc[mask_plus, "start_read"].values
    )

    # Reverse strand (-): distance = PAS position (start_read) - read genomic end (end_read)
    mask_minus = valid_reads["strand"] == "-"
    distance_df.loc[mask_minus, "distance"] = (
        valid_reads.loc[mask_minus, "start_read"].values
        - valid_reads.loc[mask_minus, "end_read"].values
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
