import logging
from typing import Optional, Tuple

import pandas as pd

from .config_manager import ConfigManager

logger = logging.getLogger(__name__)


def load_pas_sites(pas_bed_path: str) -> pd.DataFrame:
    """
    Load polyadenylation sites from a BED file.

    Parameters
    ----------
    pas_bed_path : str
        Path to polyadenylation sites BED file (chr, start, end, ..., strand format)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chr, start, end, strand
    
    Raises
    ------
    RuntimeError
        If no PAS sites are found in the file
    """
    logger.info("[1/6] Loading polyadenylation sites from BED file...")
    
    pas_df = pd.read_csv(
        pas_bed_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 5],
        names=["chr", "start", "end", "strand"],
    )
    
    if pas_df.empty:
        raise RuntimeError(f"No PAS sites found in {pas_bed_path}")
    
    logger.info(f"      ✓ Loaded {len(pas_df)} PAS sites from BED file")
    return pas_df


def get_single_pas(pas_df: pd.DataFrame) -> tuple:
    """
    Extract the single PAS site to process.
    
    Currently uses the first PAS site from the BED file.
    
    Parameters
    ----------
    pas_df : pd.DataFrame
        DataFrame with PAS sites (columns: chr, start, end, strand)
    
    Returns
    -------
    tuple
        (chr, strand, pas_position) for the PAS to process
    """
    pas_row = pas_df.iloc[0]
    pas_chr = pas_row["chr"]
    pas_strand = pas_row["strand"]
    pas_pos = pas_row["start"] if pas_strand == "-" else pas_row["end"]
    
    logger.info(f"      ✓ Processing PAS at {pas_chr}:{pas_pos} ({pas_strand})")
    
    return pas_chr, pas_strand, pas_pos
