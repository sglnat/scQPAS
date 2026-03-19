# scQPAS Code Review - Issues Found

## Summary
Found **1 critical bug** and **multiple logic issues** that will affect pipeline correctness and robustness. These range from data type errors to semantic logic problems in distance calculation.

---

## 🔴 CRITICAL BUG

### 1. **process_pas_BED.py - `get_adj_gtf()` function (Line ~160)**

**Severity**: CRITICAL - Will cause runtime failure

**Location**: [scqpas/process_pas_BED.py](scqpas/process_pas_BED.py#L160)

**Issue**: 
```python
adjusted_gtf = adjusted_gtf["start"] - 1  # BUG: Reassigns to Series instead of updating DataFrame
```

**Problem**: This line reassigns `adjusted_gtf` to a pandas Series object (the result of subtracting 1 from a column), not a DataFrame. All subsequent operations expecting a DataFrame will fail.

**Fix**: Should be:
```python
adjusted_gtf["start"] = adjusted_gtf["start"] - 1
```

**Impact**: The `get_adj_gtf()` function will crash when trying to drop columns or apply any DataFrame methods later.

---

## 🟠 LOGIC ISSUES

### 2. **calculate_distances.py - Incorrect distance calculation logic**

**Severity**: HIGH - Core pipeline metric is wrong

**Location**: [scqpas/calculate_distances.py](scqpas/calculate_distances.py#L195-L209)

**Issue**:
```python
# Calculate distance for strand '+'
mask_plus = valid_reads["strand"] == "+"
distance_df.loc[mask_plus, "distance"] = (
    valid_reads.loc[mask_plus, "end_read"].values
    - valid_reads.loc[mask_plus, "start_read"].values  # This is READ LENGTH, not distance to PAS
)

# Calculate distance for strand '-'
mask_minus = valid_reads["strand"] == "-"
distance_df.loc[mask_minus, "distance"] = (
    valid_reads.loc[mask_minus, "start_read"].values
    - valid_reads.loc[mask_minus, "end_read"].values  # This is READ LENGTH, not distance to PAS
)
```

**Problem**: 
- Current logic calculates read length (end - start), not distance from read to PAS position
- The docstring says "Computes pairwise distances between read end positions and annotated intron positions" but the code doesn't use PAS information at all
- Distance should be calculated relative to the PAS position (from `adjust_read_ends`), not just the read span
- The function receives `valid_reads` which doesn't have the adjusted coordinates from `adjust_read_ends`

**Expected Logic**:
The distance should be calculated from where the read ends to where the PAS is located. The `adjust_read_ends()` function modifies read endpoints to match PAS positions, but these adjusted coordinates are not used in the distance calculation.

**Impact**: All downstream distance metrics will be incorrect. These are the main output of the pipeline.

---

### 3. **core.py - Distance calculation uses non-adjusted read coordinates**

**Severity**: HIGH - Uses stale data for calculation

**Location**: [scqpas/core.py](scqpas/core.py#L220-L230)

**Issue**:
```python
# Step 7: CIGAR validation on rs_pas_adj (adjusted coordinates)
rs_pas_adj = adjust_read_ends(rs_pas_df)

# Step 8: But distance calculation uses valid_reads (original, non-adjusted coordinates)
valid_reads = filter_by_cigar(reads_df, rs_pas_df, reads_introns_df, cigar_df)
distances_df = calculate_distances(valid_reads, reads_introns_df)  # Uses original reads!
```

**Problem**:
1. `adjust_read_ends()` adjusts read coordinates to match PAS positions
2. These adjusted coordinates are computed but only used for CIGAR validation
3. Distance calculation receives `valid_reads` with original, non-adjusted coordinates
4. This means distances are calculated from original read positions, not from PAS-adjusted positions

**Fix**: 
- Either pass adjusted coordinates to `calculate_distances()`, or
- Rethink what "distance" should mean in this context (distance from read to nearest intron boundary? distance from read to PAS?)

**Impact**: Distances are calculated from wrong reference points.

---

### 4. **process_pas_BED.py - `assign_rs_pas()` step 7 Cartesian product issue**

**Severity**: MEDIUM - Creates incorrect read-PAS associations

**Location**: [scqpas/process_pas_BED.py](scqpas/process_pas_BED.py#L330-L340)

**Issue**:
```python
# Step 7: Expand back to read level
result_df = rs_pas_df.merge(
    reads_unique,
    on="rs_id"  # Merges on rs_id only!
)
```

**Problem**:
- After filtering for downstream PAS at the rs->transcript->PAS level, the code expands back to reads
- The merge only uses `rs_id`, which creates a Cartesian product
- If an rs_id has transcripts T1 and T2, with PAS assignments, but reads are only from T1, the merge will still assign all combinations to all reads
- This breaks the transcript-read association

**Fix**: The merge should maintain the (rs_id, transcript_id) association:
```python
result_df = rs_pas_df.merge(
    reads_unique,
    on=["rs_id", "transcript_id"]  # Add transcript_id to maintain association
)
```

**Impact**: Reads may be assigned to PAS associated with transcripts they're not mapped to.

---

### 5. **process_pas_BED.py - `adjust_read_ends()` format validation missing**

**Severity**: MEDIUM - No error handling for malformed data

**Location**: [scqpas/process_pas_BED.py](scqpas/process_pas_BED.py#L370-L380)

**Issue**:
```python
def adjust_read_ends(rs_pas_df: pd.DataFrame) -> pd.DataFrame:
    # Extract PAS components from pas_id (e.g., chr1:10465:+)
    rs_pas_df[['chr_pas', 'pos_pas', 'strand_pas']] = rs_pas_df['pas_id'].str.split(':', n=2, expand=True)
    rs_pas_df['pos_pas'] = rs_pas_df['pos_pas'].astype(int)  # Will crash if format is wrong
```

**Problem**:
- Assumes `pas_id` format is `chr:pos:strand`
- No validation of the format before splitting
- If a pas_id doesn't contain exactly 2 colons, the split will create NaN values
- The `astype(int)` conversion will fail if `pos_pas` contains non-numeric values
- No error message to help debug which row caused the problem

**Fix**: Add validation:
```python
def adjust_read_ends(rs_pas_df: pd.DataFrame) -> pd.DataFrame:
    # Validate pas_id format
    if not rs_pas_df['pas_id'].str.contains(':').all():
        raise ValueError("Invalid pas_id format. Expected 'chr:pos:strand'")
    
    parts = rs_pas_df['pas_id'].str.split(':', n=2, expand=True)
    if parts.shape[1] != 3:
        raise ValueError("pas_id must contain exactly 2 colons (format: chr:pos:strand)")
    
    rs_pas_df[['chr_pas', 'pos_pas', 'strand_pas']] = parts
    
    try:
        rs_pas_df['pos_pas'] = rs_pas_df['pos_pas'].astype(int)
    except ValueError as e:
        raise ValueError(f"pos_pas column contains non-numeric values: {e}")
```

**Impact**: Silent data corruption or cryptic error messages if PAS format is unexpected.

---

## 🟡 POTENTIAL ISSUES

### 6. **extract_reads_BAM.py - NaN handling for non-polyA reads**

**Severity**: LOW-MEDIUM - Data type inconsistency

**Location**: [scqpas/extract_reads_BAM.py](scqpas/extract_reads_BAM.py#L250-L260)

**Issue**:
```python
def get_cpa_sites(is_polyA: bool, end: int) -> Optional[int]:
    if not is_polyA:
        cpa_site = None  # Returns None
    if is_polyA:
        cpa_site = end    # Returns int
    return cpa_site
```

**Problem**:
- Returns `None` for non-polyA reads, but the column is later used in groupby operations
- In `calculate_distances.py`, NaN values in `cpa_site` may cause unexpected behavior in distance calculations
- The `dropna(subset=["cpa_site"])` in `get_cpa_sites()` function explicitly removes these, which is fine, but the inconsistency could cause issues if code changes

**Impact**: Downstream code must properly handle None/NaN values. Currently handled but fragile.

---

### 7. **extract_annotation_GTF.py - Hardcoded column indices without validation**

**Severity**: LOW-MEDIUM - Brittle to input format changes

**Location**: Multiple locations in [scqpas/extract_annotation_GTF.py](scqpas/extract_annotation_GTF.py)

**Issue**:
```python
# In extract_exons(), extract_genes(), extract_transcripts():
exons = exons[[0, 3, 4, 8, 6]].copy()  # Hardcoded indices
```

**Problem**:
- Assumes GTF file has exactly 9 columns (standard format)
- No validation that input has expected number of columns
- If GTF has extra columns or is malformed, the indexing will fail with cryptic error
- No error message indicates what went wrong

**Fix**: Add column count validation:
```python
def extract_exons(gtf_input: Union[str, pd.DataFrame], ...):
    if isinstance(gtf_input, str):
        gtf_df = load_gtf(gtf_input)
    else:
        gtf_df = gtf_input.copy()
    
    # Validate GTF format
    if gtf_df.shape[1] < 9:
        raise ValueError(f"GTF must have at least 9 columns, got {gtf_df.shape[1]}")
    
    exons = gtf_df[gtf_df[2] == "exon"]
    # ... rest of code
```

**Impact**: Poor error messages if input is malformed. Users won't know what went wrong.

---

### 8. **cli.py - No validation that BAM file is indexed**

**Severity**: LOW - Will fail at runtime with unclear error

**Location**: [scqpas/cli.py](scqpas/cli.py#L15-L35)

**Issue**:
```python
@click.option(
    "--bam",
    type=click.Path(exists=True),  # Only checks existence, not indexing
    required=True,
    help="Path to input BAM file (must be indexed)",
)
```

**Problem**:
- The docstring says BAM "must be indexed" but the validation only checks if the file exists
- If user provides an unindexed BAM, the error will occur deep in `pysam.AlignmentFile()` with a cryptic message
- No early validation to provide helpful error message

**Fix**: Add BAM index validation:
```python
import os
# Before calling run_pipeline:
if not os.path.exists(bam + ".bai"):
    raise click.ClickException(f"BAM file is not indexed. Index it with: samtools index {bam}")
```

**Impact**: Users get confusing error messages instead of helpful guidance.

---

### 9. **bedtools_intersections.py - No handling of empty outputs**

**Severity**: LOW - Edge case handling

**Location**: [scqpas/bedtools_intersections.py](scqpas/bedtools_intersections.py#L55-L70)

**Issue**:
```python
if result.stdout.strip():
    with open(output_path, "w") as f:
        f.write(result.stdout)
    return pd.read_csv(output_path, sep="\t", header=None)
else:
    return pd.DataFrame()  # Returns empty DataFrame with no columns
```

**Problem**:
- When bedtools returns no intersections, an empty DataFrame is returned with no columns
- Downstream code expects columns [0, 1, 2, ...] for column assignment
- If no intersections exist, column assignment like `df[[0, 1, 2, ...]].copy()` in `core.py` will crash

**Example of crash**:
```python
# In core.py Step 5:
reads_transcripts_df = run_bedtools_intersect(...)  # Returns DataFrame with 0 rows, 0 columns
reads_transcripts_df = reads_transcripts_df[[0, 1, 2, 3, 4, 5, 9, 12]].copy()  # KeyError: 0
```

**Fix**: Return DataFrame with expected column structure:
```python
if result.stdout.strip():
    # ... existing code
else:
    # Return empty DataFrame with expected columns
    return pd.DataFrame()
```

And in `core.py`, handle empty results properly:
```python
if not reads_transcripts_df.empty:
    # ... process
else:
    logger.warning("No intersections found")
    distances_df = pd.DataFrame()
```

**Impact**: Pipeline crashes if no read-transcript intersections exist, rather than producing empty output.

---

## 📋 Summary Table

| Issue | File | Severity | Type |
|-------|------|----------|------|
| DataFrame reassignment to Series | process_pas_BED.py | 🔴 CRITICAL | Bug |
| Distance calculation uses read length not PAS distance | calculate_distances.py | 🟠 HIGH | Logic |
| Distance calculation uses non-adjusted reads | core.py | 🟠 HIGH | Logic |
| Cartesian product in read expansion | process_pas_BED.py | 🟠 MEDIUM | Logic |
| No format validation for pas_id | process_pas_BED.py | 🟠 MEDIUM | Robustness |
| NaN handling in CPA sites | extract_reads_BAM.py | 🟡 LOW | Edge case |
| Hardcoded column indices | extract_annotation_GTF.py | 🟡 LOW | Robustness |
| No BAM index validation | cli.py | 🟡 LOW | Usability |
| Empty bedtools output handling | bedtools_intersections.py | 🟡 LOW | Edge case |

---

## 🔧 Recommended Fix Priority

1. **FIRST**: Fix the critical DataFrame reassignment bug (Issue #1)
2. **SECOND**: Fix the distance calculation logic (Issues #2, #3)
3. **THIRD**: Fix the Cartesian product issue (Issue #4)
4. **FOURTH**: Add data validation and error handling (Issues #5-9)

