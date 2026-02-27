"""
Bedtools intersection utilities for scQPAS.

This module handles running bedtools intersect with flexible flags and managing
BED file inputs/outputs.
"""

import logging
import os
import subprocess
import pandas as pd

logger = logging.getLogger(__name__)


def run_bedtools_intersect(
    bed_a_df, bed_b_df, tmpdir, name_a="a", name_b="b", flags=None
):
    """
    Run bedtools intersect with flexible flags.

    Parameters
    ----------
    bed_a_df : pd.DataFrame
        BED format DataFrame (chr, start, end, ...)
    bed_b_df : pd.DataFrame
        BED format DataFrame
    tmpdir : str
        Temporary directory for intermediate files
    name_a, name_b : str
        Names for temporary files
    flags : list, optional
        bedtools flags (e.g., ['-wa', '-s', '-f', '1.0'])
        Default: ['-wo'] (write overlapping positions from both)

    Returns
    -------
    pd.DataFrame
        Intersection results with numeric column names (0, 1, 2, ...)

    Raises
    ------
    FileNotFoundError
        If bedtools is not found
    RuntimeError
        If intersection fails
    """

    bed_a_path = os.path.join(tmpdir, f"bed_{name_a}.bed")
    bed_b_path = os.path.join(tmpdir, f"bed_{name_b}.bed")
    output_path = os.path.join(tmpdir, f"intersection_{name_a}_{name_b}.bed")

    # Write temporary BED files
    bed_a_df.to_csv(bed_a_path, sep="\t", header=False, index=False)
    bed_b_df.to_csv(bed_b_path, sep="\t", header=False, index=False)

    # Run bedtools intersect with custom flags
    try:

        cmd = ["bedtools", "intersect"] + flags + ["-a", bed_a_path, "-b", bed_b_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if result.stdout.strip():
            with open(output_path, "w") as f:
                f.write(result.stdout)
            return pd.read_csv(output_path, sep="\t", header=None)
        else:
            return pd.DataFrame()
    except FileNotFoundError:
        raise FileNotFoundError(
            "bedtools not found. Please install bedtools: conda install -c bioconda bedtools"
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bedtools error: {e.stderr}")
