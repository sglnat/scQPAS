import pandas as pd
import logging
from typing import Optional
from .config_manager import ConfigManager
import csv
import gc

logger = logging.getLogger(__name__)


def add_read_set_ids(
    reads_transcripts_df: pd.DataFrame, reads_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Add read set IDs (rs_id) to reads_transcripts_df by tracing read_ids back to original reads.

    Creates a composed ID (rs_id) from cell barcode and UMI in the format: CB_UMI.
    Merges this information back onto the reads_transcripts_df based on read_id.

    Parameters
    ----------
    reads_transcripts_df : pd.DataFrame
        DataFrame from bedtools intersect with columns: chr, start, end, read_id, dummy, strand, transcript_id
    reads_df : pd.DataFrame
        Original reads DataFrame from extract_reads with columns including: read_id, UMI, CB

    Returns
    -------
    pd.DataFrame
        reads_transcripts_df with additional 'rs_id' column (format: CB_UMI)
    """
    # Extract only read_id, UMI, and CB from original reads
    read_set_info = reads_df[["read_id", "UMI", "CB"]].copy()

    # Create rs_id as CB_UMI
    read_set_info["rs_id"] = read_set_info["CB"] + "_" + read_set_info["UMI"]

    # Keep only read_id and rs_id, then merge on read_id onto reads_transcripts_df
    read_set_info = read_set_info[["read_id", "rs_id"]]
    reads_transcripts_df = reads_transcripts_df.merge(
        read_set_info, on="read_id", how="left"
    )

    return reads_transcripts_df


def filter_by_transcript_association(
    reads_transcripts_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Filter reads based on transcript_id association within read sets.

    For each read set (rs_id) and transcript_id combination, checks if ALL reads in that
    read set are assigned to this transcript. Only keeps (rs_id, transcript_id) pairs where
    every read in the read set has that transcript assignment.

    Parameters
    ----------
    reads_transcripts_df : pd.DataFrame
        DataFrame with columns including: rs_id, read_id, transcript_id
        (each row represents one read-transcript assignment; reads can have multiple assignments)

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame keeping only reads with transcript assignments valid for their read set
    """
    # For each (rs_id, transcript_id) pair, count how many unique reads have this assignment
    reads_per_pair = reads_transcripts_df.groupby(["rs_id", "transcript_id"])[
        "read_id"
    ].transform("nunique")

    # For each rs_id, count total unique reads in that read set
    total_reads_per_rs = reads_transcripts_df.groupby("rs_id")["read_id"].transform(
        "nunique"
    )

    # Keep only rows where ALL reads in the read set are assigned to this transcript
    reads_transcripts_df = reads_transcripts_df[reads_per_pair == total_reads_per_rs]

    return reads_transcripts_df.reset_index(drop=True)


def _adjust_terminal_exons(
    gtf: pd.DataFrame, extension_length: int = 1000
) -> pd.DataFrame:
    """
    Extend terminal exons to capture extended 3' UTR regions.

    For each transcript, modifies the coordinates of the terminal (last) exon
    and associated transcript/gene records. Extension direction depends on strand:
    - Forward strand (+): extend end coordinate by extension_length bp (rightmost exon)
    - Reverse strand (-): extend start coordinate by -extension_length bp (leftmost exon)

    Parameters
    ----------
    gtf : pd.DataFrame
        GTF DataFrame with columns: seqname, source, feature, start, end, strand,
        transcript_id, gene_id
    extension_length : int, optional
        Length in base pairs to extend terminal exons. Default: 1000

    Returns
    -------
    pd.DataFrame
        Modified GTF DataFrame with extended terminal exons (preserves row order)
    """
    gtf = gtf.copy()

    # Identify terminal exons: forward strand by max end, reverse by min start
    # Work only with exon rows from GTF
    exon_df = gtf[gtf["feature"] == "exon"]

    # Forward strand: get indices of terminal exons (max end per transcript)
    fw_exons = exon_df[exon_df["strand"] == "+"]
    fw_terminal_indices = fw_exons.loc[
        fw_exons.groupby("transcript_id")["end"].idxmax()
    ].index
    fw_terminal_exons = fw_exons.loc[fw_terminal_indices]

    # Reverse strand: get indices of terminal exons (min start per transcript)
    rv_exons = exon_df[exon_df["strand"] == "-"]
    rv_terminal_indices = rv_exons.loc[
        rv_exons.groupby("transcript_id")["start"].idxmin()
    ].index
    rv_terminal_exons = rv_exons.loc[rv_terminal_indices]

    # Get transcript_ids and gene_ids of terminal exons
    fw_transcript_ids = fw_terminal_exons["transcript_id"].unique()
    rv_transcript_ids = rv_terminal_exons["transcript_id"].unique()
    fw_gene_ids = fw_terminal_exons["gene_id"].unique()
    rv_gene_ids = rv_terminal_exons["gene_id"].unique()

    # Update terminal exons using their actual indices
    gtf.loc[fw_terminal_indices, "end"] += extension_length
    gtf.loc[rv_terminal_indices, "start"] -= extension_length

    # Update associated transcript records
    gtf.loc[
        (gtf["feature"] == "transcript")
        & (gtf["transcript_id"].isin(fw_transcript_ids)),
        "end",
    ] += extension_length
    gtf.loc[
        (gtf["feature"] == "transcript")
        & (gtf["transcript_id"].isin(rv_transcript_ids)),
        "start",
    ] -= extension_length

    # Update associated gene records
    gtf.loc[
        (gtf["feature"] == "gene") & (gtf["gene_id"].isin(fw_gene_ids)), "end"
    ] += extension_length
    gtf.loc[
        (gtf["feature"] == "gene") & (gtf["gene_id"].isin(rv_gene_ids)), "start"
    ] -= extension_length

    return gtf


def get_adj_gtf(
    input_file,
    output_file: Optional[str] = None,
    extension_length: Optional[int] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Load GTF file or DataFrame, extend terminal exons, and optionally save adjusted GTF.

    Reads a GTF file or uses a GTF DataFrame, parses attributes, adjusts terminal exons,
    and optionally writes the modified GTF to file.

    Parameters
    ----------
    input_file : str or pd.DataFrame
        Path to input GTF file OR GTF DataFrame. If string, reads from file. If DataFrame, uses directly.
    output_file : str, optional
        Path to output GTF file with extended terminal exons. If None, returns DataFrame only (in-memory).
    extension_length : int, optional
        Length in base pairs to extend terminal exons. If None, retrieves from config_manager. Default from config: 1000
    config_manager : ConfigManager, optional
        Configuration manager for accessing default extension_length. Default: None

    Returns
    -------
    pd.DataFrame
        The adjusted GTF DataFrame
    """
    # Get extension_length from config if not provided
    if extension_length is None:
        if config_manager is None:
            config_manager = ConfigManager()
        extension_length = config_manager.get("annotation", "terminal_exon_extension")

    # Load GTF from file or use provided DataFrame
    if isinstance(input_file, str):
        gtf = pd.read_csv(input_file, sep="\t", comment="#", header=None)
    else:
        gtf = input_file

    # Set columns if not already set (from file read)
    if len(gtf.columns) == 9:
        gtf.columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

    # Extract transcript_id and gene_id for grouping
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    gtf["gene_id"] = gtf["attribute"].str.extract(r'gene_id "([^"]+)"')

    adjusted_gtf = _adjust_terminal_exons(gtf, extension_length=extension_length)

    # Drop temporary columns
    adjusted_gtf = adjusted_gtf.drop(columns=["transcript_id", "gene_id"])

    if output_file is not None:
        adjusted_gtf.to_csv(
            output_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
        )

    return adjusted_gtf


def load_pas(
    pas_file_path: str, stringency: int, region: Optional[str] = None
) -> pd.DataFrame:
    """
    Load a PAS BED file into memory with optional coordinate filtering.

    Parameters
    ----------
    pas_file_path : str
        Path to PAS BED file
    stringency : int
        Minimum stringency level for PAS. Only PAS with stringency >= this value will be retained.
    region : str, optional
        Optional genomic region for filtering in format "chr:start-end" (e.g., "chr12:15671-26783")
        or just chromosome (e.g., "chr12"). If None, no coordinate filtering is applied.
        Default: None (loads entire file)

    Returns
    -------
    pd.DataFrame
        PAS DataFrame with columns: chr, start, end, pas_id, gex, strand
        Filtered by stringency and optionally by region coordinates.
    """

    pas = pd.read_csv(pas_file_path, sep="\t", header=None)
    pas.columns = [
        "chr",
        "start",
        "end",
        "pas_id",
        "gex",
        "strand",
        "tissues",
        "protocols",
        "stringency",
        "gen_class",
        "polyA_signal",
    ]
    pas = pas[pas["stringency"] >= stringency]

    # Optional region-based filtering
    if region is not None:
        # Parse region string
        if ":" in region:
            # Format: "chr12:15671-26783"
            chr_part, coord_part = region.split(":")
            start_coord, end_coord = coord_part.split("-")
            start_coord = int(start_coord)
            end_coord = int(end_coord)
        else:
            # Format: "chr12" (chromosome only)
            chr_part = region
            start_coord = None
            end_coord = None

        # Filter by chromosome
        pas = pas[pas["chr"] == chr_part]

        # Filter by coordinates if provided
        if start_coord is not None and end_coord is not None:
            pas = pas[(pas["start"] >= start_coord) & (pas["end"] <= end_coord)]

    return pas[["chr", "start", "end", "pas_id", "gex", "strand"]]


def _get_terminal_reads_per_rs(reads_transcripts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify terminal reads for each read set (rs_id).

    For each rs_id, identifies the terminal (most 3') read based on strand:
    - Forward strand (+): rightmost read (maximum end_read coordinate)
    - Reverse strand (-): leftmost read (minimum start_read coordinate)

    Parameters
    ----------
    reads_transcripts_df : pd.DataFrame
        DataFrame with columns must include: rs_id, chr_read, start_read, end_read, strand.
        Must be pre-filtered to unique read_ids per rs_id.

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per rs_id containing terminal read information:
        rs_id, chr_read, terminal_start, terminal_end, strand
    """
    df = reads_transcripts_df

    # Split into forward and reverse strand groups
    fw_mask = df["strand"] == "+"
    rv_mask = df["strand"] == "-"

    # Forward strand: group by rs_id and get row with maximum end_read
    fw_df = df[fw_mask].loc[df[fw_mask].groupby("rs_id")["end_read"].idxmax()].copy()
    fw_df.rename(
        columns={"start_read": "terminal_start", "end_read": "terminal_end"},
        inplace=True,
    )
    fw_df = fw_df[["rs_id", "chr_read", "terminal_start", "terminal_end", "strand"]]

    # Reverse strand: group by rs_id and get row with minimum start_read
    rv_df = df[rv_mask].loc[df[rv_mask].groupby("rs_id")["start_read"].idxmin()].copy()
    rv_df.rename(
        columns={"start_read": "terminal_start", "end_read": "terminal_end"},
        inplace=True,
    )
    rv_df = rv_df[["rs_id", "chr_read", "terminal_start", "terminal_end", "strand"]]

    # Concatenate both strand groups
    return pd.concat([fw_df, rv_df], ignore_index=False)


def assign_rs_pas(
    reads_transcripts_df: pd.DataFrame,
    pas_transcript_df: pd.DataFrame,
    config_manager: Optional[ConfigManager] = None,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Assign PAS to read sets based on transcript association and terminal read position.

    Merge-based approach working at rs/transcript level, then expanding to read level:
    1. Get terminal read position per rs_id (vectorized, one row per rs)
    2. Get unique rs->transcript associations
    3. Merge rs->transcript with transcript->pas (creates rs->transcript->pas grid)
    4. Filter PAS downstream of terminal read
    5. Expand back to read level using read_id mapping

    A PAS is considered "downstream" if:
    - Forward strand (+): PAS end_pas >= terminal read end_read (3' boundary inclusive)
    - Reverse strand (-): PAS start_pas <= terminal read start_read (3' boundary inclusive)

    Parameters
    ----------
    reads_transcripts_df : pd.DataFrame
        DataFrame with columns: chr_read, start_read, end_read, read_id,
        strand, transcript_id, rs_id
    pas_transcript_df : pd.DataFrame
        DataFrame with columns: chr_pas, start_pas, end_pas, pas_id, gex,
        strand_pas, transcript_id, gene_id
    config_manager : ConfigManager, optional
        Configuration manager for accessing bedtools_score setting. Default: None
    output_csv : str, optional
        Path to output CSV file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        One row per read-PAS association with columns:
        rs_id, read_id, transcript_id, pas_id, and coordinates from both inputs
    """
    reads_df = reads_transcripts_df
    # cols: "chr_read", "start_read", "end_read", "read_id", "dummy", "strand", "transcript_id", "gene_id", "rs_id"
    pas_df = pas_transcript_df
    # cols: "chr_pas", "start_pas", "end_pas", "pas_id", "gex", "strand_pas", "transcript_id", "gene_id"

    # Step 1: Get unique read_id per rs_id before computing terminal positions
    # This serves as both the source for terminal reads and the mapping for later expansion
    reads_unique = reads_df[
        ["rs_id", "read_id", "chr_read", "start_read", "end_read", "strand"]
    ].drop_duplicates(subset=["rs_id", "read_id"])

    # Step 2: Get terminal read info per rs_id
    terminal_reads_df = _get_terminal_reads_per_rs(reads_unique)

    # Step 3: Create rs->transcript associations
    rs_transcript_df = reads_df[["rs_id", "transcript_id"]].drop_duplicates()

    # Step 4: Merge rs->transcript with transcript->pas
    # This creates a Cartesian product: for each (rs_id, transcript_id) pair,
    # we get all PAS associated with that transcript_id
    # Example: if rs1 has tr1, and tr1 has pas1, pas2 -> result has (rs1, tr1, pas1) and (rs1, tr1, pas2)
    rs_pas_df = rs_transcript_df.merge(
        pas_df[
            ["transcript_id", "pas_id", "chr_pas", "start_pas", "end_pas", "strand_pas"]
        ],
        on="transcript_id",
    )

    # Step 5: Merge with terminal reads to add position filtering info
    rs_pas_df = rs_pas_df.merge(terminal_reads_df, on="rs_id")

    # Step 6: Filter PAS that are downstream of terminal read
    # Forward strand: PAS end >= terminal read end; Reverse strand: PAS start <= terminal read start
    downstream = (rs_pas_df["chr_pas"] == rs_pas_df["chr_read"]) & (
        (
            (rs_pas_df["strand"] == "+")
            & (rs_pas_df["end_pas"] >= rs_pas_df["terminal_end"])
        )
        | (
            (rs_pas_df["strand"] == "-")
            & (rs_pas_df["start_pas"] <= rs_pas_df["terminal_start"])
        )
    )

    rs_pas_df = rs_pas_df[downstream][["rs_id", "transcript_id", "pas_id"]]

    # Step 7: Expand back to read level using mapping from reads_unique
    # For each (rs_id, transcript_id, pas_id) at rs level, assign to all reads in that rs
    result_df = rs_pas_df.merge(reads_unique, on="rs_id")

    # Step 8: Reorder columns in bed order, add dummy col
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")

    result_df["dummy"] = bedtools_score
    result_df = result_df[
        [
            "chr_read",
            "start_read",
            "end_read",
            "read_id",
            "dummy",
            "strand",
            "transcript_id",
            "pas_id",
            "rs_id",
        ]
    ]

    result_df = result_df.reset_index(drop=True)

    if output_csv is not None:
        result_df.to_csv(output_csv, index=False)

    return result_df


def adjust_read_ends(
    rs_pas_df: pd.DataFrame, output_csv: Optional[str] = None
) -> pd.DataFrame:
    """
    Adjust read end coordinates to reflect PAS position for terminal reads.

    For each read associated with a PAS, updates the read's end coordinate to match the PAS position.

    Parameters
    ----------
    rs_pas_df : pd.DataFrame
        DataFrame with columns: rs_id, read_id, transcript_id, pas_id, chr_read, start_read, end_read, strand
    output_csv : str, optional
        Path to output CSV file. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        DataFrame with adjusted read coordinates for terminal reads associated with PAS
    """

    # Extract PAS components from pas_id (e.g., chr1:10465:+)
    rs_pas_df[["chr_pas", "pos_pas", "strand_pas"]] = rs_pas_df["pas_id"].str.split(
        ":", n=2, expand=True
    )
    rs_pas_df["pos_pas"] = rs_pas_df["pos_pas"].astype(int)

    # Create masks for strands
    fw_mask = rs_pas_df["strand"] == "+"
    rv_mask = rs_pas_df["strand"] == "-"

    # Forward strand: set end_read to PAS position
    rs_pas_df.loc[fw_mask, "end_read"] = rs_pas_df.loc[fw_mask, "pos_pas"]

    # Reverse strand: set start_read to PAS position - 1
    rs_pas_df.loc[rv_mask, "start_read"] = rs_pas_df.loc[rv_mask, "pos_pas"] - 1

    if output_csv is not None:
        rs_pas_df.to_csv(output_csv, index=False)

    return rs_pas_df


def filter_polyA_by_pas(
    polyA_reads_df: pd.DataFrame,
    pas_df: pd.DataFrame,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Filter polyA read sets by PAS and assign pas_ids to detected CPA sites.

    **Algorithm:**
    1. Add rs_id (CB_UMI) to each polyA read using add_read_set_ids()
    2. For each RS, identify all unique cpa_sites from polyA reads in that RS
    3. For each cpa_site, check if it overlaps with any atlas PAS entry
    4. Only keep (RS, cpa_site) combinations that match atlas PAS
    5. Expand reads so each read has a row for every valid cpa_site in its RS
       (allows reads with different cpa_sites within same RS to calculate distances to all)
    6. Assign matching pas_id in format "chr:pos:strand" from atlas PAS

    **Key insight:** A single RS may have multiple cpa_sites (if polyA reads detected
    different cleavage positions). Each read must be evaluated against all valid
    cpa_sites in its RS, producing up to reads_per_RS * unique_cpa_sites rows.

    Parameters
    ----------
    polyA_reads_df : pd.DataFrame
        polyA reads from extract_reads(). Columns must include:
        read_id, UMI, CB, chr, start, end, strand, cpa_site
    pas_df : pd.DataFrame
        Pre-filtered atlas PAS. Columns must include:
        chr, start, end, pas_id, strand
        (typically from load_pas() with region and stringency filtering already applied)
    output_csv : str, optional
        Path to save output CSV. If None, returns DataFrame only (in-memory).

    Returns
    -------
    pd.DataFrame
        Indexed polyA reads with one row per (read, cpa_site) pair.
        Columns: read_id, rs_id, chr, start, end, strand, cpa_site, pas_id, UMI, CB

        Only includes reads where cpa_site matches an atlas PAS on the same chromosome
        and strand.
    """

    # Step 1: Add rs_id (CB_UMI) directly to input dataframe
    polyA_reads_df["rs_id"] = polyA_reads_df["CB"] + "_" + polyA_reads_df["UMI"]

    # Step 2: Validate chr/strand consistency and filter to valid RS's
    rs_consistency = polyA_reads_df.groupby("rs_id").agg(
        {"chr": "nunique", "strand": "nunique"}
    )
    valid_rs_ids = rs_consistency[
        (rs_consistency["chr"] == 1) & (rs_consistency["strand"] == 1)
    ].index
    del rs_consistency

    polyA_reads = polyA_reads_df[polyA_reads_df["rs_id"].isin(valid_rs_ids)]

    # Get unique cpa_sites per RS (excluding NaN), along with chr and strand
    rs_cpa_sites = (
        polyA_reads.groupby("rs_id")
        .agg(
            {
                "cpa_site": lambda x: x.dropna().unique(),
                "chr": "first",
                "strand": "first",
            }
        )
        .reset_index()
    )

    # Step 3: Explode cpa_sites array into individual rows
    cpa_df = rs_cpa_sites.explode("cpa_sites").reset_index(drop=True)
    cpa_df.rename(columns={"cpa_sites": "cpa_site"}, inplace=True)
    del rs_cpa_sites

    # Step 4: Match cpa_sites to overlapping PAS
    # Merge on chr/strand to get all PAS candidates, then filter to overlapping only
    rs_cpa_pas_df = cpa_df.merge(
        pas_df[["chr", "strand", "start", "end", "pas_id"]],
        on=["chr", "strand"],
        how="left",
    )

    # Keep overlapping PAS (or NaN if no PAS exists for this chr/strand)
    rs_cpa_pas_df = rs_cpa_pas_df[
        (rs_cpa_pas_df["start"].isna())
        | (
            (rs_cpa_pas_df["start"] <= rs_cpa_pas_df["cpa_site"])
            & (rs_cpa_pas_df["cpa_site"] <= rs_cpa_pas_df["end"])
        )
    ][["rs_id", "cpa_site", "pas_id"]].drop_duplicates()
    del cpa_df

    # Step 5: Expand to read level using merge
    result_df = polyA_reads[
        ["read_id", "rs_id", "chr", "start", "end", "strand", "UMI", "CB"]
    ].merge(rs_cpa_pas_df, on="rs_id")
    del polyA_reads, rs_cpa_pas_df

    # Reorder columns
    result_df = result_df[
        [
            "read_id",
            "rs_id",
            "chr",
            "start",
            "end",
            "strand",
            "cpa_site",
            "pas_id",
            "UMI",
            "CB",
        ]
    ]

    if result_df.empty:
        logger.warning("No polyA RS had cpa_sites matching atlas PAS")
    else:
        logger.info(
            f"Filtered polyA RS to {len(result_df['rs_id'].unique())} RS with "
            f"{len(result_df['cpa_site'].unique())} unique valid CPA sites"
        )

    if output_csv is not None:
        result_df.to_csv(output_csv, index=False)

    gc.collect()
    return result_df
