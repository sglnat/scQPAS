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
    Count occurrences of each cleavage/polyadenylation site by chromosome.

    Groups reads by chromosome and CPA site position to identify all putative
    cleavage sites and their supporting read counts.

    Parameters
    ----------
    reads_df : pd.DataFrame
        DataFrame containing read information with columns 'chr' and 'cpa_site'.
        Typically the output from extract_reads() in extract_reads_BAM.py.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['chr', 'cpa_site', 'supporting_read_count']
        showing the number of reads supporting each (chr, cpa_site) pair.
        Sorted by supporting_read_count in descending order.
    """
    # Remove rows with missing CPA sites
    df_valid = reads_df.dropna(subset=["cpa_site"])

    if df_valid.empty:
        raise RuntimeError("No valid CPA sites found in reads")

    # Count occurrences of each (chr, cpa_site) pair
    cpa_counts = (
        df_valid.groupby(["chr", "cpa_site"])
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

    Selects the (chr, cpa_site) pair with the highest supporting read count.
    This represents the most probable cleavage site in the dataset.

    Parameters
    ----------
    cpa_sites_df : pd.DataFrame
        DataFrame with columns ['chr', 'cpa_site', 'supporting_read_count']
        as returned by get_cpa_sites(). Should already be sorted by
        supporting_read_count in descending order.

    Returns
    -------
    tuple or None
        (chr, cpa_site_position) for the CPA site with highest read count,
        or None if no CPA sites are found
    """
    if cpa_sites_df.empty:
        return None

    # Get the first row (highest supporting read count due to sorting)
    best_row = cpa_sites_df.iloc[0]
    best_chr = best_row["chr"]
    best_pos = best_row["cpa_site"]

    return (best_chr, best_pos)


def detect_best_cpa_from_reads(
    reads_df: pd.DataFrame,
) -> tuple:
    """
    Automatically detect the best CPA site from read data.

    This function analyzes the reads to identify the cleavage/polyadenylation site
    with the most supporting evidence (reads with detected polyA tails).

    Parameters
    ----------
    reads_df : pd.DataFrame
        DataFrame from extract_reads() containing columns: read_id, chr, cpa_site,
        is_polyA, and other read information.

    Returns
    -------
    tuple
        (chr, pas_pos, strand) where:
        - chr: chromosome of the detected CPA site
        - pas_pos: genomic position of the detected CPA site
        - strand: strand information ('.' if not determinable from data,
                  or '+'/'-' if available from reads at this position)

    Raises
    ------
    RuntimeError
        If no polyA reads are found or if no valid CPA sites can be detected
    """
    # Filter for polyA reads only (these are the ones with detected CPA sites)
    polyA_reads = reads_df[reads_df["is_polyA"] == True].copy()

    if polyA_reads.empty:
        raise RuntimeError(
            "No polyA reads found in BAM file. "
            "Cannot detect CPA site from data. "
            "Check if your BED file PAS coordinates are correct and overlap with reads."
        )

    # Get CPA site distribution
    cpa_sites = get_cpa_sites(polyA_reads)

    # Get the best (most supported) CPA
    best_chr, best_cpa_pos = best_cpa(cpa_sites)

    # Determine strand from reads at this CPA position
    reads_at_cpa = polyA_reads[
        (polyA_reads["chr"] == best_chr) & (polyA_reads["cpa_site"] == best_cpa_pos)
    ]

    # Get the most common strand at this position
    strand_counts = reads_at_cpa["strand"].value_counts()
    best_strand = strand_counts.index[0] if not strand_counts.empty else "."

    logger.info(
        f"Detected best CPA site: {best_chr}:{best_cpa_pos} ({best_strand}) "
        f"with {reads_at_cpa.shape[0]} polyA-containing reads"
    )

    return best_chr, best_cpa_pos, best_strand


def calculate_distances(
    valid_reads: pd.DataFrame, intronic_reads: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate distances from reads to cleavage/polyadenylation sites.

    Computes pairwise distances between read end positions and annotated intron
    positions to quantify the distance between where the read ends and where
    cleavage is predicted to occur.

    Parameters
    ----------
    valid_reads : pd.DataFrame
        DataFrame with validated reads including columns:
        'read_id', 'transcript_id', 'intronic', 'CIGAR_N'
    intronic_reads : pd.DataFrame
        Intersection results with columns:
        'chr_read', 'start_read', 'end_read', 'read_id', 'length_intron', etc.

    Returns
    -------
    pd.DataFrame
        DataFrame with distance calculations and supporting information
    """

    intronic_reads["transcript_id"] = intronic_reads["intron_id"].str.split("_").str[0]
    # print(intronic_reads)

    # Create a new DataFrame for storing distances
    distance_df = pd.DataFrame()
    distance_df[["read_id", "transcript_id", "intronic", "CIGAR_N"]] = valid_reads[
        ["read_id", "transcript_id", "intronic", "CIGAR_N"]
    ]
    # print(distance_df)

    # Merging to include total intron lengths in the distance_df
    distance_df = distance_df.merge(
        intronic_reads.groupby(["read_id", "transcript_id"])["length_intron"]
        .sum()
        .reset_index(),
        on=["read_id", "transcript_id"],
        how="left",
    )
    distance_df.rename(columns={"length_intron": "total_intron_len"}, inplace=True)

    # Calculate distance based on strand information
    distance_df["distance"] = None

    # Calculate distance for strand '+'
    mask_plus = valid_reads["strand"] == "+"
    distance_df.loc[mask_plus, "distance"] = (
        valid_reads.loc[mask_plus, "end"].values
        - valid_reads.loc[mask_plus, "start"].values
    )

    # Calculate distance for strand '-'
    mask_minus = valid_reads["strand"] == "-"
    distance_df.loc[mask_minus, "distance"] = (
        valid_reads.loc[mask_minus, "start"].values
        - valid_reads.loc[mask_minus, "end"].values
    )

    # Subtract total_intron_len from distance for intronic reads
    intronic_mask = valid_reads["intronic"] == True
    distance_df.loc[intronic_mask, "distance"] -= distance_df.loc[
        intronic_mask, "total_intron_len"
    ]

    # print(distance_df)

    return distance_df
