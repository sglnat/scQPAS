import logging
from typing import Optional, Union

import pandas as pd
import numpy as np

from .config_manager import ConfigManager

logger = logging.getLogger(__name__)


def load_gtf(gtf_file_path: str) -> pd.DataFrame:
    """
    Load a GTF file into memory with standard column names.

    Parameters
    ----------
    gtf_file_path : str
        Path to GTF file

    Returns
    -------
    pd.DataFrame
        GTF DataFrame with named columns: seqname, source, feature, start, end, score, strand, frame, attribute
    """
    gtf = pd.read_csv(gtf_file_path, sep="\t", comment="#", header=None)
    gtf.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    return gtf


def extract_exons(
    gtf_input: Union[str, pd.DataFrame],
    exons_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Extract exons from a GTF file or DataFrame.

    Parameters
    ----------
    gtf_input : str or pd.DataFrame
        Path to GTF file or GTF DataFrame. If string, reads from file. If DataFrame, uses directly.
    exons_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager (not currently used by this function)

    Returns
    -------
    pd.DataFrame
        Exons with columns: chr, start, end, transcript_id, strand, gene_id
    """
    # Load GTF if path is provided, otherwise use DataFrame directly
    if isinstance(gtf_input, str):
        gtf_df = load_gtf(gtf_input)
    else:
        gtf_df = gtf_input
    
    # Ensure named columns if not already set
    if not hasattr(gtf_df, 'columns') or len(gtf_df.columns) == 0 or gtf_df.columns[0] == 0:
        gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    
    exons = gtf_df[gtf_df["feature"] == "exon"][["seqname", "start", "end", "attribute", "strand"]]
    exons.columns = ["chr", "start", "end", "attributes", "strand"]
    exons["start"] = exons["start"] - 1  # Convert start to 0-based
    exons["transcript_id"] = exons["attributes"].str.extract(
        r'transcript_id\s+"([^"]+)"'
    )
    exons["gene_id"] = exons["attributes"].str.extract(r'gene_id\s+"([^"]+)"')
    exons.drop(columns="attributes", inplace=True)

    # Get bedtools score from config
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")

    exons["dummy"] = bedtools_score
    # Final reordering - single operation to avoid intermediate DFs
    exons = exons[["chr", "start", "end", "transcript_id", "dummy", "strand", "gene_id"]]

    if exons_output is not None:
        exons.to_csv(exons_output, sep="\t", header=False, index=False)

    return exons


def calculate_introns(
    exons: pd.DataFrame,
    introns_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Calculate introns from exons.

    Returns DataFrame with transcript-specific introns, defined as the regions
    between exons of the same transcript.

    Parameters
    ----------
    exons : pd.DataFrame
        Exons DataFrame from extract_exons()
    introns_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager (not currently used by this function)

    Returns
    -------
    pd.DataFrame
        Introns with columns: chr, start, end, intron_id, length_intron, strand, gene_id
    """

    grouped_exons = exons.groupby("transcript_id")

    introns_list = []

    for transcript_id, group in grouped_exons:
        group = group.sort_values(by="start")
        exons_starts = group["start"].to_list()
        exons_ends = group["end"].to_list()

        # Calculate intron positions
        for i in range(len(exons_starts) - 1):
            intron_start = exons_ends[i]
            intron_end = exons_starts[i + 1]
            intron_length = intron_end - intron_start
            introns_list.append(
                {
                    "chr": group["chr"].iloc[0],
                    "start": intron_start,
                    "end": intron_end,
                    "intron_id": f"{transcript_id}_{i}",
                    "length_intron": intron_length,
                    "strand": group["strand"].iloc[0],
                    "gene_id": group["gene_id"].iloc[0],
                }
            )

    introns_df = pd.DataFrame(introns_list)

    if introns_output is not None:
        introns_df.to_csv(introns_output, sep="\t", header=False, index=False)

    return introns_df


def extract_genes(
    gtf_input: Union[str, pd.DataFrame],
    genes_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Extract genes from a GTF file or DataFrame.

    Parameters
    ----------
    gtf_input : str or pd.DataFrame
        Path to GTF file or GTF DataFrame. If string, reads from file. If DataFrame, uses directly.
    genes_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager (not currently used by this function)

    Returns
    -------
    pd.DataFrame
        Genes with columns: chr, start, end, gene_id, strand
    """
    # Load GTF if path is provided, otherwise use DataFrame directly
    if isinstance(gtf_input, str):
        gtf_df = load_gtf(gtf_input)
    else:
        gtf_df = gtf_input
    
    # Ensure named columns if not already set
    if not hasattr(gtf_df, 'columns') or len(gtf_df.columns) == 0 or gtf_df.columns[0] == 0:
        gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    
    genes = gtf_df[gtf_df["feature"] == "gene"][["seqname", "start", "end", "attribute", "strand"]]
    genes.columns = ["chr", "start", "end", "attributes", "strand"]
    genes["start"] = genes["start"] - 1  # Convert start to 0-based
    genes["gene_id"] = genes["attributes"].str.extract(r'gene_id\s+"([^"]+)"')
    genes.drop(columns="attributes", inplace=True)

    # Get bedtools score from config
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")

    genes["dummy"] = bedtools_score
    # Final reordering - single operation to avoid intermediate DFs
    genes = genes[["chr", "start", "end", "gene_id", "dummy", "strand"]]

    if genes_output is not None:
        genes.to_csv(genes_output, sep="\t", header=False, index=False)

    return genes


def extract_transcripts(
    gtf_input: Union[str, pd.DataFrame],
    transcripts_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Extract transcripts from a GTF file or DataFrame.

    Parameters
    ----------
    gtf_input : str or pd.DataFrame
        Path to GTF file or GTF DataFrame. If string, reads from file. If DataFrame, uses directly.
    transcripts_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager (not currently used by this function)

    Returns
    -------
    pd.DataFrame
        Transcripts with columns: chr, start, end, transcript_id, strand
    """
    # Load GTF if path is provided, otherwise use DataFrame directly
    if isinstance(gtf_input, str):
        gtf_df = load_gtf(gtf_input)
    else:
        gtf_df = gtf_input
    
    # Ensure named columns if not already set
    if not hasattr(gtf_df, 'columns') or len(gtf_df.columns) == 0 or gtf_df.columns[0] == 0:
        gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    
    transcripts = gtf_df[gtf_df["feature"] == "transcript"][["seqname", "start", "end", "attribute", "strand"]]
    transcripts.columns = ["chr", "start", "end", "attributes", "strand"]
    transcripts["start"] = transcripts["start"] - 1  # Convert start to 0-based
    transcripts["transcript_id"] = transcripts["attributes"].str.extract(r'transcript_id\s+"([^"]+)"')
    transcripts["gene_id"] = transcripts["attributes"].str.extract(r'gene_id\s+"([^"]+)"')
    transcripts.drop(columns="attributes", inplace=True)

    # Get bedtools score from config
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")

    transcripts["dummy"] = bedtools_score
    # Final reordering - single operation to avoid intermediate DFs
    transcripts = transcripts[["chr", "start", "end", "transcript_id", "dummy", "strand", "gene_id"]]

    if transcripts_output is not None:
        transcripts.to_csv(transcripts_output, sep="\t", header=False, index=False)

    return transcripts


def reads_to_bed(
    df_input: pd.DataFrame,
    bed_output: Optional[str] = None,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Convert reads DataFrame to BED format.

    Parameters
    ----------
    df_input : pd.DataFrame
        Reads DataFrame with columns: chr, start, end, read_id, strand
    bed_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager for accessing bedtools_score setting. Default: None

    Returns
    -------
    pd.DataFrame
        BED format DataFrame with columns: chr, start, end, read_id, dummy, strand
    """
    df = df_input.copy()
    
    # Convert start to 0-based coordinates for BED format
    df["start"] = df["start"] - 1
    
    # Get bedtools score from config
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")
    
    # Add dummy score column
    df["dummy"] = bedtools_score
    
    # Reorder columns to BED format
    df = df[["chr", "start", "end", "read_id", "dummy", "strand"]]
    
    if bed_output is not None:
        df.to_csv(bed_output, sep="\t", header=False, index=False)
    
    return df