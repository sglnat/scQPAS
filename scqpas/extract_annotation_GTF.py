import logging
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def extract_exons(gtf_file_path, exons_output=None, config_manager=None):
    """
    Extract exons from a GTF file.

    Parameters
    ----------
    gtf_file_path : str
        Path to GTF file
    exons_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager (not currently used by this function)

    Returns
    -------
    pd.DataFrame
        Exons with columns: chr, start, end, transcript_id, strand, gene_id
    """

    gtf_df = pd.read_csv(gtf_file_path, sep="\t", comment="#", header=None)
    exons = gtf_df[gtf_df[2] == "exon"]
    exons = exons[
        [0, 3, 4, 8, 6]
    ].copy()  # 0: chromosome, 3: start, 4: end, 6: strand, 8: attributes
    exons.columns = ["chr", "start", "end", "attributes", "strand"]
    exons["start"] = exons["start"] - 1  # Convert start to 0-based
    exons["transcript_id"] = exons["attributes"].str.extract(
        r'transcript_id\s+"([^"]+)"'
    )
    exons["gene_id"] = exons["attributes"].str.extract(r'gene_id\s+"([^"]+)"')
    exons.drop(columns="attributes", inplace=True)
    exons = exons[["chr", "start", "end", "transcript_id", "strand", "gene_id"]]

    if exons_output is not None:
        exons.to_csv(
            exons_output,
            sep="\t",
            header=False,
            index=False,
            columns=["chr", "start", "end", "transcript_id", "strand", "gene_id"],
        )

    return exons


def calculate_introns(exons, introns_output=None, config_manager=None):
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


def get_bed_from_df(df_input, chr, cpa_site, bed_output=None, config_manager=None):
    """
    Create BED format for reads with PAS-specific coordinate adjustment.

    Parameters
    ----------
    df_input : pd.DataFrame or str
        Reads DataFrame or path to CSV file with columns: chr, start, end, strand, read_id
    chr : str
        Target chromosome (e.g., 'chr12')
    cpa_site : int
        Cleavage/polyadenylation site position (genomic coordinate)
    bed_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager for accessing bedtools_score setting. Default: None

    Returns
    -------
    pd.DataFrame
        BED format DataFrame with columns: chr, start, end, read_id, dummy, strand

    Notes
    -----
    For strand-specific coordinate adjustment:
    - Forward strand: Use region from read start to PAS position
    - Reverse strand: Use region from PAS position to read end
    The dummy column is populated with bedtools_score from config (default: 1000)
    """
    # Handle both DataFrame and file path inputs
    if isinstance(df_input, str):
        df = pd.read_csv(df_input, sep=",")
    else:
        df = df_input.copy()

    bed = []

    for idx, row in df.iterrows():
        if row["strand"] == "+" and row["chr"] == chr:
            start = row["start"]
            end = cpa_site + 1  # BED end is exclusive
        elif row["strand"] == "-" and row["chr"] == chr:
            start = cpa_site  # BED start is inclusive
            end = row["end"]
        else:
            continue  # Skip rows with an undefined strand

        if start < end:
            bed.append([row["chr"], start, end, row["read_id"], row["strand"]])

    bed_df = pd.DataFrame(bed, columns=["chr", "start", "end", "read_id", "strand"])

    # Get bedtools score from config with fallback to default
    if config_manager:
        bedtools_score = config_manager.get("output", "bedtools_score", 1000)
    else:
        bedtools_score = 1000

    bed_df.insert(4, "dummy", bedtools_score)

    if bed_output is not None:
        bed_df.to_csv(bed_output, sep="\t", header=False, index=False)

    return bed_df


def extract_genes(gtf_file_path, genes_output=None, config_manager=None):
    """
    Extract genes from a GTF file.

    Parameters
    ----------
    gtf_file_path : str
        Path to GTF file
    genes_output : str, optional
        Path to output BED file. If None, returns DataFrame only (in-memory).
    config_manager : ConfigManager, optional
        Configuration manager for accessing bedtools_score setting. Default: None

    Returns
    -------
    pd.DataFrame
        Genes with columns: chr, start, end, gene_id, dummy, strand
        The dummy column is populated with bedtools_score from config (default: 1000)
    """

    gtf_df = pd.read_csv(gtf_file_path, sep="\t", comment="#", header=None)
    genes = gtf_df[gtf_df[2] == "gene"]
    genes = genes[[0, 3, 4, 8, 6]].copy()
    genes.columns = ["chr", "start", "end", "attributes", "strand"]
    genes["start"] = genes["start"] - 1
    genes["gene_id"] = genes["attributes"].str.extract(r'gene_id\s+"([^"]+)"')
    genes.drop(columns="attributes", inplace=True)
    genes = genes[["chr", "start", "end", "gene_id", "strand"]]

    # Get bedtools score from config with fallback to default
    if config_manager:
        bedtools_score = config_manager.get("output", "bedtools_score", 1000)
    else:
        bedtools_score = 1000

    genes.insert(4, "dummy", bedtools_score)

    if genes_output is not None:
        genes.to_csv(
            genes_output,
            sep="\t",
            header=False,
            index=False,
            columns=["chr", "start", "end", "gene_id", "dummy", "strand"],
        )

    return genes
