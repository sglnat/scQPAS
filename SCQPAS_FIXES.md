# scQPAS - Detailed Fixes for Critical Issues

## Issue #1: DataFrame Reassignment Bug (CRITICAL)

### Current Code
**File**: `scqpas/process_pas_BED.py` (~line 160)

```python
def get_adj_gtf(input_file: str, output_file: Optional[str] = None) -> pd.DataFrame:
    # ... code ...
    
    adjusted_gtf = adjusted_gtf["start"] - 1  # BUG: Returns Series, not DataFrame!
    
    # Drop temporary columns
    adjusted_gtf = adjusted_gtf.drop(columns=["transcript_id", "gene_id"])  # CRASH HERE
```

### Fixed Code
```python
def get_adj_gtf(input_file: str, output_file: Optional[str] = None) -> pd.DataFrame:
    gtf = pd.read_csv(input_file, sep="\t", comment="#", header=None)
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

    adjusted_gtf = adjust_terminal_exons(gtf)

    # FIX: Use assignment operator to update the column, not reassign the whole DataFrame
    adjusted_gtf["start"] = adjusted_gtf["start"] - 1  # Convert to 0-based coordinates for BED format

    # Drop temporary columns
    adjusted_gtf = adjusted_gtf.drop(columns=["transcript_id", "gene_id"])

    if output_file is not None:
        adjusted_gtf.to_csv(
            output_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
        )

    return adjusted_gtf
```

### What Changed
- **Line changed**: `adjusted_gtf = adjusted_gtf["start"] - 1` → `adjusted_gtf["start"] = adjusted_gtf["start"] - 1`
- **Impact**: The variable now remains a DataFrame instead of being converted to a Series

---

## Issue #2 & #3: Distance Calculation Logic (HIGH)

### Current Code
**File**: `scqpas/calculate_distances.py` (~lines 175-209)

```python
def calculate_distances(
    valid_reads: pd.DataFrame, intronic_reads: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate distances from reads to cleavage/polyadenylation sites.
    """
    distance_df = valid_reads[["read_id", "rs_id", "transcript_id", "pas_id", "intronic", "CIGAR_N"]].copy()

    distance_df = distance_df.merge(
        intronic_reads.groupby(["read_id", "transcript_id", "pas_id"])["length_intron"]
        .sum()
        .reset_index(),
        on=["read_id", "transcript_id", "pas_id"],
        how="left",
    )
    distance_df.rename(columns={"length_intron": "total_intron_len"}, inplace=True)

    distance_df["distance"] = None

    # BUG: This calculates READ LENGTH, not distance to PAS
    mask_plus = valid_reads["strand"] == "+"
    distance_df.loc[mask_plus, "distance"] = (
        valid_reads.loc[mask_plus, "end_read"].values
        - valid_reads.loc[mask_plus, "start_read"].values  # This is read span!
    )

    mask_minus = valid_reads["strand"] == "-"
    distance_df.loc[mask_minus, "distance"] = (
        valid_reads.loc[mask_minus, "start_read"].values
        - valid_reads.loc[mask_minus, "end_read"].values
    )

    intronic_mask = valid_reads["intronic"] == True
    distance_df.loc[intronic_mask, "distance"] -= distance_df.loc[
        intronic_mask, "total_intron_len"
    ]

    return distance_df
```

### Problem Analysis

The current logic calculates read length (end - start coordinates), not the distance from a read to a cleavage site. Looking at the pipeline:

1. `adjust_read_ends()` is called on `rs_pas_df` to set read end coordinates to match PAS position
2. But `adjust_read_ends()` result is NOT passed to `calculate_distances()`
3. Instead, `valid_reads` (original, non-adjusted) is passed
4. So distances are calculated from original read boundaries, not PAS positions

### Questions to Clarify

Before implementing the fix, the team needs to clarify:
- **What should "distance" actually measure?**
  - Option A: Distance from read 3' boundary to PAS position?
  - Option B: Something else based on intron information?
  
- **Why is `adjust_read_ends()` called if its results aren't used?**
  - Should adjusted coordinates be used in distance calculation?
  - Or is there a different purpose for `adjust_read_ends()`?

### Proposed Fix (Option A: Use adjusted coordinates)

Modify `core.py` Step 8 to pass adjusted coordinates:

```python
# Step 7: BEDTOOLS INTERSECTIONS: READS-INTRONS
logger.info("[7/8] Intersecting reads with annotated introns...")

rs_pas_adj = adjust_read_ends(rs_pas_df)

# ... bedtools intersection code ...

# Step 8: CALCULATE DISTANCES
logger.info("[8/8] Calculating distances from CPA sites...")

valid_reads = filter_by_cigar(
    reads_df, rs_pas_df, reads_introns_df, cigar_df
)

# IMPORTANT: Merge valid_reads with adjusted coordinates
valid_reads_adj = valid_reads.merge(
    rs_pas_adj[["read_id", "start_read", "end_read"]],
    on="read_id",
    how="left",
    suffixes=("_orig", "_adj")
)

# Use adjusted read ends for distance calculation
distances_df = calculate_distances(valid_reads_adj, reads_introns_df)
```

And modify `calculate_distances.py`:

```python
def calculate_distances(
    valid_reads: pd.DataFrame, intronic_reads: pd.DataFrame, use_adjusted: bool = True
) -> pd.DataFrame:
    """
    Calculate distances from reads to cleavage/polyadenylation sites.
    
    Parameters
    ----------
    valid_reads : pd.DataFrame
        DataFrame with validated reads. If use_adjusted=True, should contain
        '_adj' suffixed columns for PAS-adjusted coordinates.
    ...
    """
    distance_df = valid_reads[["read_id", "rs_id", "transcript_id", "pas_id", "intronic", "CIGAR_N"]].copy()

    # Merge intron lengths
    distance_df = distance_df.merge(
        intronic_reads.groupby(["read_id", "transcript_id", "pas_id"])["length_intron"]
        .sum()
        .reset_index(),
        on=["read_id", "transcript_id", "pas_id"],
        how="left",
    )
    distance_df.rename(columns={"length_intron": "total_intron_len"}, inplace=True)

    # Choose which coordinates to use
    if use_adjusted and "end_read_adj" in valid_reads.columns:
        start_col = "start_read_adj"
        end_col = "end_read_adj"
    else:
        start_col = "start_read"
        end_col = "end_read"

    distance_df["distance"] = 0

    # For forward strand: distance from adjusted read end to PAS region
    # PAS was already set as read end by adjust_read_ends, so distance = 0 for terminal reads
    # But we should calculate distance from original read to PAS position
    mask_plus = valid_reads["strand"] == "+"
    distance_df.loc[mask_plus, "distance"] = (
        valid_reads.loc[mask_plus, end_col].values
        - valid_reads.loc[mask_plus, start_col].values
    )

    # For reverse strand
    mask_minus = valid_reads["strand"] == "-"
    distance_df.loc[mask_minus, "distance"] = (
        valid_reads.loc[mask_minus, start_col].values
        - valid_reads.loc[mask_minus, end_col].values
    )

    # Subtract intron lengths for intronic reads
    intronic_mask = valid_reads["intronic"] == True
    distance_df.loc[intronic_mask, "distance"] -= distance_df.loc[
        intronic_mask, "total_intron_len"
    ]

    return distance_df
```

**Question**: The logic still seems unclear. More context is needed about what "distance" actually represents in the context of polyA site detection. Is it:
- Distance from 3' end of read to PAS?
- Read length minus intron length?
- Something else?

---

## Issue #4: Cartesian Product in Read Expansion (HIGH)

### Current Code
**File**: `scqpas/process_pas_BED.py` (~lines 330-340)

```python
def assign_rs_pas(...):
    # ... steps 1-6 ...
    
    # Step 7: Expand back to read level using mapping from reads_unique
    result_df = rs_pas_df.merge(
        reads_unique,
        on="rs_id"  # BUG: Only merges on rs_id, creates Cartesian product
    )
```

### Problem Example

```
rs_pas_df:
  rs_id  transcript_id  pas_id
  rs1    T1             pas1
  rs1    T1             pas2

reads_unique:
  rs_id  read_id  transcript_id
  rs1    r1       T1
  rs1    r2       T1

Result of on="rs_id" merge (WRONG):
  rs_id  transcript_id_x  pas_id  read_id  transcript_id_y
  rs1    T1               pas1    r1       T1
  rs1    T1               pas1    r2       T1      <- Wrong! r2 gets pas1
  rs1    T1               pas2    r1       T1
  rs1    T1               pas2    r2       T1      <- Wrong! r2 gets pas2

Result of on=["rs_id", "transcript_id"] merge (CORRECT):
  rs_id  transcript_id  pas_id  read_id
  rs1    T1             pas1    r1
  rs1    T1             pas1    r2
  rs1    T1             pas2    r1
  rs1    T1             pas2    r2
```

### Fixed Code

```python
def assign_rs_pas(
    reads_transcripts_df: pd.DataFrame, 
    pas_transcript_df: pd.DataFrame,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
    """
    Assign PAS to read sets based on transcript association and terminal read position.
    """
    reads_df = reads_transcripts_df.copy()
    pas_df = pas_transcript_df.copy()
    
    # Step 1: Get unique read_id per rs_id
    reads_unique = reads_df[["rs_id", "read_id", "chr_read", "start_read", "end_read", "strand"]].drop_duplicates(subset=["rs_id", "read_id"])
    
    # Step 2: Get terminal read info per rs_id
    terminal_reads_df = get_terminal_reads_per_rs(reads_unique)
    
    # Step 3: Create rs->transcript associations
    rs_transcript_df = reads_df[["rs_id", "transcript_id"]].drop_duplicates()
    
    # Step 4: Merge rs->transcript with transcript->pas
    rs_pas_df = rs_transcript_df.merge(
        pas_df[["transcript_id", "pas_id", "chr_pas", "start_pas", "end_pas", "strand_pas"]],
        on="transcript_id"
    )
    
    # Step 5: Merge with terminal reads
    rs_pas_df = rs_pas_df.merge(terminal_reads_df, on="rs_id")
    
    # Step 6: Filter PAS downstream of terminal read
    downstream = (rs_pas_df["chr_pas"] == rs_pas_df["chr_read"]) & (
        ((rs_pas_df["strand"] == "+") & (rs_pas_df["end_pas"] >= rs_pas_df["terminal_end"])) |
        ((rs_pas_df["strand"] == "-") & (rs_pas_df["start_pas"] <= rs_pas_df["terminal_start"]))
    )
    rs_pas_df = rs_pas_df[downstream][["rs_id", "transcript_id", "pas_id"]]   
    
    # Step 7: Expand back to read level - FIX: Merge on both rs_id AND transcript_id
    result_df = rs_pas_df.merge(
        reads_unique[["rs_id", "read_id", "chr_read", "start_read", "end_read", "strand"]],
        on=["rs_id", "transcript_id"]  # FIX: Added transcript_id to maintain association
    )
    
    # Step 8: Add metadata and reorder columns
    if config_manager is None:
        config_manager = ConfigManager()
    bedtools_score = config_manager.get("output", "bedtools_score")

    result_df["dummy"] = bedtools_score
    result_df = result_df[["chr_read", "start_read", "end_read", "read_id", "dummy", "strand", "transcript_id", "pas_id", "rs_id"]]
        
    return result_df.reset_index(drop=True)
```

### What Changed
- **Line changed**: `on="rs_id"` → `on=["rs_id", "transcript_id"]`
- **Impact**: Maintains the correct transcript association for each read, avoiding cross-assignment between transcripts

---

## Issue #5: Format Validation for pas_id (MEDIUM)

### Current Code
**File**: `scqpas/process_pas_BED.py` (~lines 370-380)

```python
def adjust_read_ends(rs_pas_df: pd.DataFrame) -> pd.DataFrame:
    rs_pas_df[['chr_pas', 'pos_pas', 'strand_pas']] = rs_pas_df['pas_id'].str.split(':', n=2, expand=True)
    rs_pas_df['pos_pas'] = rs_pas_df['pos_pas'].astype(int)  # Silent failure if format wrong
```

### Fixed Code

```python
def adjust_read_ends(rs_pas_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adjust read end coordinates to reflect PAS position for terminal reads.

    Parameters
    ----------
    rs_pas_df : pd.DataFrame
        DataFrame with columns: rs_id, read_id, transcript_id, pas_id, chr_read, start_read, end_read, strand
        Note: pas_id must have format 'chr:pos:strand'

    Returns
    -------
    pd.DataFrame
        DataFrame with adjusted read coordinates for terminal reads associated with PAS
    """
    
    # Validate pas_id format before processing
    invalid_pas = rs_pas_df[~rs_pas_df['pas_id'].str.contains(':')]
    if not invalid_pas.empty:
        logger.error(f"Found {len(invalid_pas)} pas_id entries without ':' separator")
        logger.error(f"Example invalid pas_id: {invalid_pas['pas_id'].iloc[0]}")
        raise ValueError("pas_id must have format 'chr:pos:strand' with ':' separators")
    
    # Extract PAS components from pas_id (e.g., chr1:10465:+)
    parts = rs_pas_df['pas_id'].str.split(':', n=2, expand=True)
    
    if parts.shape[1] != 3:
        invalid_count = len(rs_pas_df[~rs_pas_df['pas_id'].str.contains(':', regex=False)])
        raise ValueError(
            f"pas_id format error: expected 'chr:pos:strand' format. "
            f"Found {invalid_count} entries that don't follow this format. "
            f"Example: {rs_pas_df['pas_id'].iloc[0]}"
        )
    
    rs_pas_df[['chr_pas', 'pos_pas', 'strand_pas']] = parts
    
    # Validate that pos_pas can be converted to integer
    try:
        rs_pas_df['pos_pas'] = rs_pas_df['pos_pas'].astype(int)
    except ValueError as e:
        invalid_rows = rs_pas_df[~rs_pas_df['pos_pas'].str.isnumeric()]
        logger.error(f"Found {len(invalid_rows)} non-numeric position values")
        logger.error(f"Example invalid position: {invalid_rows['pos_pas'].iloc[0]}")
        raise ValueError(
            f"pas_id position values must be numeric. "
            f"Invalid values found. Check pas_id format in BED file."
        ) from e
    
    # Create masks for strands
    fw_mask = rs_pas_df["strand"] == "+"
    rv_mask = rs_pas_df["strand"] == "-"
    
    # Forward strand: set end_read to PAS position
    rs_pas_df.loc[fw_mask, "end_read"] = rs_pas_df.loc[fw_mask, "pos_pas"]
    
    # Reverse strand: set start_read to PAS position - 1
    rs_pas_df.loc[rv_mask, "start_read"] = rs_pas_df.loc[rv_mask, "pos_pas"] - 1
    
    return rs_pas_df
```

### What Changed
- Added validation before splitting
- Provides clear error messages with examples
- Helps users debug if PAS format is incorrect

