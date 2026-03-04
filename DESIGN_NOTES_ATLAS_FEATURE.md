# Design Notes: Genome-Wide PAS Atlas Feature

**Status:** Planning phase (not yet implemented)

**Related Code:** 
- [scqpas/calculate_distances.py](scqpas/calculate_distances.py) - `detect_best_cpa_from_reads()`, `get_cpa_sites()`, `best_cpa()`
- [scqpas/core.py](scqpas/core.py) - Step 2: CPA site detection
- [scqpas/cli.py](scqpas/cli.py) - `--pas` argument (currently deprecated but reserved for atlas)

---

## 1. Current State (V0: Single CPA Site)

### Single PAS Workflow
```
BAM file → Extract reads (detect polyA) → Identify best CPA from reads 
  → GTF annotation → Calculate distances to single PAS → Output CSV
```

### Key Properties
- **One distance per read:** Each read has a single distance value to the detected CPA site
- **Data-driven PAS:** No user-specified PAS; detected automatically from polyA read enrichment
- **Strand-aware:** CPA position accounts for transcript strand (start vs end)
- **No filtering:** All polyA reads included in CPA site frequency analysis

### Output Structure (Current)
```csv
read_id,transcript_id,intronic,CIGAR_N,total_intron_len,distance,pas_coordinate
READ001,ENST0001,True,2,450,1250,42000000
READ002,ENST0001,False,0,0,625,42000000
...
```

---

## 2. Future State (V1: Multiple PAS per Read)

### Extended Workflow with Atlas
```
┌─────────────────┐
│  BAM file       │
└────────┬────────┘
         │ Extract reads + detect polyA
         ▼
┌─────────────────────────────────┐
│  Reads DataFrame                │
│  (chr, pos, strand, is_polyA)   │
└────────┬────────────────────────┘
         │ (Step 2a: Identify best CPA from reads)
         │          [for reporting/validation]
         ▼
┌──────────────────────────────────┐
│ ANNOTATION                       │
├──────────────────────────────────┤
│ Step 3: GTF annotation           │
│ Step 4: CIGAR-derived introns    │
└────────┬─────────────────────────┘
         │
         ▼
    ┌────────────────────┐
    │  Load PAS ATLAS    │  ◄── NEW: Genome-wide BED file
    │  (--pas-atlas)     │        with all known PAS
    │  e.g., Atlas_hg38  │
    │  Millions of sites │
    └────────┬───────────┘
             │ Filter PAS candidates
             │ for each read
             ▼
    ┌──────────────────────────────────┐
    │ PAS Filtering Strategy:           │
    │ • By chromosome (match read)      │
    │ • By strand (match read strand)   │
    │ • By distance from read 3' end    │
    │   (e.g., within ±5kb)             │
    └────────┬─────────────────────────┘
             │
             ▼
  ┌──────────────────────────────────┐
  │ Reads-to-PAS Mapping             │
  │ (Many-to-many: N reads × M PAS)  │
  │ filtered_pas_per_read = [        │
  │   (pos1, dist1),                 │
  │   (pos2, dist2),                 │
  │   (pos3, dist3)                  │
  │ ]                                │
  └────────┬─────────────────────────┘
           │ Bedtools intersections
           │ (same as before)
           ▼
  ┌────────────────────────────────────┐
  │ Calculate Distances (per read,     │
  │ per candidate PAS)                 │
  │ [replicate for M candidate PAS]    │
  └────────┬───────────────────────────┘
           │
           ▼
  ┌─────────────────────────────────────┐
  │ Output: Multiple rows per read      │
  │ (one per candidate PAS)             │
  └─────────────────────────────────────┘
```

### Key Changes from V0 to V1

#### 1. **Input Data Structure**
- **Old:** Optional `--pas` BED file (single/few sites, deprecated)
- **New:** `--pas-atlas` BED file (genome-wide, millions of sites)
- **Format:** Standard BED format
  ```
  chr1  1000000  1000100  PAS001  0  +
  chr1  1005000  1005100  PAS002  0  +
  ...
  ```

#### 2. **PAS Filtering Strategy**

The most critical design decision: **How to select candidate PAS for each read?**

##### Option A: Distance-Based Filtering (RECOMMENDED)
```python
def filter_candidate_pas(read, pas_atlas_df, 
                         distance_threshold=5000):
    """
    For each read, select all PAS:
    1. On same chromosome
    2. On same strand
    3. Within distance_threshold bp downstream (3' end)
    """
    candidates = pas_atlas_df[
        (pas_atlas_df['chr'] == read['chr']) &
        (pas_atlas_df['strand'] == read['strand']) &
        (pas_atlas_df['start'] >= read['end']) &  # downstream
        (pas_atlas_df['start'] <= read['end'] + distance_threshold)
    ]
    return candidates
```

**Advantages:**
- Biologically sensible (PAS should be near 3' end)
- Computationally efficient (reduces candidate set)
- Parameterizable (adjust distance threshold)

**Considerations:**
- Distance threshold depends on gene structure, splicing patterns
- May need different thresholds for + vs - strand
- Should validate empirically on test dataset

##### Option B: Gene-Based Filtering (Alternative)
```python
def filter_candidate_pas_by_gene(read, genes_df, pas_atlas_df):
    """
    For each read:
    1. Find intersecting gene(s)
    2. Select all PAS for that gene
    3. Filter by strand
    """
    # Would require gene annotation linking to PAS
    # (possible if PAS atlas provides gene_id)
```

**Advantages:**
- More specific to gene structure
- Avoids spurious distant PAS

**Challenges:**
- Requires PAS-to-gene mapping in atlas
- Reads may map to multiple genes
- More complex logic

---

#### 3. **Data Structure Changes in Core Pipeline**

##### A. Reads DataFrame (Current, Unchanged)
```python
# From extract_reads_BAM.py
reads_df.columns = [
    'read_id', 'UMI', 'CB', 'chr', 'start', 'end', 'strand',
    'CIGAR', 'is_polyA', 'len_pA', 'cpa_site', ... other cols
]
```

**New in V1:** 
- Add column: `candidate_pas` = list of filtered PAS coordinates
  ```python
  reads_df['candidate_pas'] = reads_df.apply(
      lambda row: filter_candidate_pas(row, pas_atlas_df),
      axis=1
  )
  ```
- OR: Create separate `reads_to_pas_mapping` DataFrame
  ```python
  reads_to_pas_mapping = pd.DataFrame({
      'read_id': [read_id, read_id, read_id, ...],
      'pas_coordinate': [pos1, pos2, pos3, ...],
      'pas_index': [idx1, idx2, idx3, ...]  # for tracing in atlas
  })
  ```

##### B. Pipeline Flow Modification

**Current Step 2 (Single CPA):**
```python
# core.py, Step 2
pas_chr, pas_pos, pas_strand = detect_best_cpa_from_reads(reads_df)
```

**V1 Step 2a (Best CPA - for reference/reporting):**
```python
# Keep this for validation and reporting
pas_chr, pas_pos, pas_strand = detect_best_cpa_from_reads(reads_df)
logger.info(f"Most common CPA in dataset: {pas_chr}:{pas_pos}")
```

**V1 Step 2b (Load Atlas & Filter Candidates):**
```python
# NEW: Load genome-wide PAS atlas
pas_atlas_df = load_pas_atlas(pas_atlas_path)

# NEW: Filter candidates for each read
reads_to_pas_mapping = filter_pas_candidates_per_read(
    reads_df, pas_atlas_df, 
    distance_threshold=5000
)
logger.info(f"Generated {len(reads_to_pas_mapping)} "
            f"read-to-PAS candidate pairs")

# NEW: Validate that we have candidates
if reads_to_pas_mapping.empty:
    raise RuntimeError(
        "No candidate PAS found for any reads. "
        "Check: chromosome naming, strand, distance threshold."
    )
```

##### C. Bedtools Intersections (Modified)

**Current (V0):**
```python
# Step 5a: Convert single read set to BED
reads_bed_df = get_bed_from_df(reads_df, pas_pos, chr=pas_chr)

# Step 5b: Single intersection
reads_genes_df = run_bedtools_intersect(reads_bed_df, genes_bed_df, ...)
```

**V1 (Multiple PAS):**
```python
# Step 5: Intersect reads with genes (UNCHANGED - per read, not per PAS)
reads_genes_df = run_bedtools_intersect(reads_bed_df, genes_bed_df, ...)

# Step 6: For each candidate PAS, calculate distances
# (distances already per-read, now iterating over candidate PAS)

for idx, (read_id, pas_row) in reads_to_pas_mapping.iterrows():
    pas_pos = pas_row['pas_coordinate']
    
    # Calculate distance: same as V0 distance calculation
    # but now inside a loop/apply over all candidate PAS
    ...
```

---

#### 4. **Distance Calculation (Modified)**

**Current `calculate_distances()` signature:**
```python
def calculate_distances(
    valid_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate distances from reads to single PAS.
    
    Returns:
        distance_df with one row per read
    """
```

**V1 Proposed Changes:**

Option A: **Extend existing function (modify signature)**
```python
def calculate_distances(
    valid_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame,
    pas_coordinates: Union[int, list, pd.Series] = None
) -> pd.DataFrame:
    """
    Calculate distances from reads to PAS (single or multiple).
    
    Parameters:
    - valid_reads: DataFrame of reads (per-read basis)
    - intronic_reads: bedtools intersection results  
    - pas_coordinates: 
        * int: single PAS (current behavior)
        * list/Series: multiple PAS (returns replicated rows)
    
    Returns:
        distance_df with:
        - V0: 1 row per read
        - V1: M rows per read (one per candidate PAS)
    """
```

Option B: **New function for atlas mode + wrapper**
```python
def calculate_distances_multi_pas(
    valid_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame,
    reads_to_pas_mapping: pd.DataFrame  # new: specifies PAS for each read
) -> pd.DataFrame:
    """
    Calculate distances from reads to multiple PAS.
    """
    
def calculate_distances(
    valid_reads: pd.DataFrame,
    intronic_reads: pd.DataFrame,
    pas_coordinates: Union[int, pd.DataFrame]  # supports both modes
) -> pd.DataFrame:
    """Wrapper that dispatches to single or multi-PAS calculator."""
```

**Recommendation:** Option B (separate function) for clarity and maintainability.

---

#### 5. **Output Format Changes**

**Current (V0):**
```csv
read_id,transcript_id,intronic,CIGAR_N,total_intron_len,distance,pas_coordinate
READ001,ENST0001,True,2,450,1250,42000000
READ002,ENST0001,False,0,0,625,42000000
READ003,ENST0002,True,2,300,850,42000000
```

**V1 Output (Multiple rows per read):**
```csv
read_id,transcript_id,intronic,CIGAR_N,total_intron_len,distance,pas_coordinate,pas_index
READ001,ENST0001,True,2,450,1250,42000000,PAS001
READ001,ENST0001,True,2,450,1200,42005000,PAS002
READ001,ENST0001,True,2,450,1100,42010000,PAS003
READ002,ENST0001,False,0,0,625,42000000,PAS001
READ002,ENST0001,False,0,0,575,42005000,PAS002
READ003,ENST0002,True,2,300,850,43000000,PAS101
```

**New Columns:**
- `pas_index`: Identifier from atlas (enables tracing back to detailed PAS info)
- `pas_id`: Optional human-readable name (e.g., "chr22:42000001-42000100")

**Additional Output File (Optional):**
```
distances_per_read_per_pas.csv  # current format
distances_summary.csv            # summary stats per read:
                                 #   min_distance, max_distance, #_candidates
```

---

## 3. Implementation Roadmap

### Phase 1: Preparation (Current - Done ✓)
- [x] Refactor `get_cpa_sites()` to include chr information
- [x] Implement `detect_best_cpa_from_reads()` for V0 single-PAS mode
- [x] Make `--pas` argument optional (reserve for V1)
- [x] Document design in this file

### Phase 2: Infrastructure (Before Atlas Implementation)
- [ ] Add `load_pas_atlas()` function in `process_pasfile.py`
- [ ] Implement `filter_pas_candidates_per_read()` with configurable distance threshold
- [ ] Create `reads_to_pas_mapping` construction logic
- [ ] Add unit tests for PAS filtering on known test data

### Phase 3: Core Implementation (Main Atlas Feature)
- [ ] Implement `calculate_distances_multi_pas()` in `calculate_distances.py`
- [ ] Modify pipeline step 2 to optionally load atlas
- [ ] Update bedtools intersection logic (if needed)
- [ ] Implement output CSV with multiple rows per read
- [ ] Add per-read summary statistics

### Phase 4: Configuration & CLI
- [ ] Add `--pas-atlas` CLI option (genome-wide atlas file)
- [ ] Add `--pas-distance-threshold` option (e.g., 5000 bp)  
- [ ] Update help text and examples
- [ ] Add mode selection: `--mode single|atlas` (or auto-detect)

### Phase 5: Testing & Validation
- [ ] Unit tests for PAS filtering logic
- [ ] Integration tests with real atlas file
- [ ] Benchmark performance (atlas loading, filtering, distance calc)
- [ ] Validation against known PAS sites

### Phase 6: Documentation
- [ ] Update README with atlas usage examples
- [ ] Create tutorial notebook with atlas workflow
- [ ] Document distance-threshold tuning strategy

---

## 4. Technical Considerations

### A. Performance & Scalability

**Current (V0):** 
- 1 bedtools intersection → 1 output CSV
- Fast (minutes for typical 1M reads)

**V1 Atlas Mode:**
- Load millions of PAS from atlas (RAM usage)
- Per-read filtering (potentially M × N operations)
- M bedtools intersections (one per candidate PAS)? OR single intersection with replicated rows?

**Optimization Strategy:**
```python
# Recommended: Single bedtools intersection, replicate rows in pandas

# Instead of:
for pas in candidate_pas_list:
    result = bedtools_intersect(reads, ...with pas...)  # M separate calls

# Do:
result = bedtools_intersect(reads, genes, ...)  # 1 call
# Then replicate output rows for each candidate PAS
distances = []
for pas in candidate_pas_list:
    distances.append(result.copy())
    distances[-1]['pas_coordinate'] = pas
distances_df = pd.concat(distances)
```

**Memory Implications:**
- Atlas loading: 10-50 MB (typical hg38 atlas with 1M+ PAS)
- Reads-to-PAS mapping: 5-100 MB (depending on candidate filtering)
- Output CSV: 100 MB - 1 GB (if M candidates × N reads × multiplier)

**Mitigation:**
- Implement chunked processing for large M values
- Consider stratified output (one CSV per chromosome)
- Option to write to database or parquet format instead of CSV

### B. Chromosome & Coordinate System

**Critical:** Ensure consistency in:
1. Chromosome naming (chr22 vs 22, chrMT vs MT)
2. Coordinate system (0-based vs 1-based)
3. Reference genome version (hg38, hg19, etc.)

**Recommendation:**
- Validate chromosome names match between BAM, GTF, and atlas
- Add explicit check at pipeline start:
  ```python
  assert_chromosome_consistency(bam_path, gtf_path, atlas_path)
  ```

### C. Strand Handling

**Current (V0):**
- Single PAS per read
- Strand already incorporated in distance calculation

**V1 Atlas:**
- Filter candidate PAS by strand match (read strand == PAS strand)
- Example: 3' end of + strand read should use PAS on + strand
- Verify logic in `filter_pas_candidates_per_read()`

---

## 5. Future Enhancements (Post-V1)

### A. Score-Based PAS Selection
- If atlas provides experimental support score (e.g., PolyASite Atlas scores)
- Instead of distance threshold, use score threshold
- Select top N PAS by score for each read

```python
def filter_by_score(read, pas_atlas_df, score_threshold=0.5):
    candidates = pas_atlas_df[
        (pas_atlas_df['chr'] == read['chr']) &
        (pas_atlas_df['score'] >= score_threshold)
    ]
    return candidates
```

### B. Gene-Specific PAS (if atlas provides gene annotation)
- Link reads to genes via bedtools intersection
- Select PAS within detected gene boundaries
- More conservative than distance-based

### C. Sample/Tissue-Specific Atlas
- Different tissues have different dominant PAS
- Allow atlas filtering by tissue type
- Prioritize tissue-matched PAS

### D. Machine Learning Integration
- Train model to predict PAS probability from read features
- Use ML scores instead of raw count/distance
- Incorporate polyA tail length, GC content, etc.

---

## 6. Testing Strategy

### Unit Tests Needed
```python
# test_calculate_distances.py

def test_filter_candidate_pas_by_distance():
    """Test PAS filtering with various distance thresholds"""
    
def test_calculate_distances_multi_pas_replicates_rows():
    """Verify output has M × N rows for M candidates × N reads"""
    
def test_chromosome_consistency_check():
    """Verify error raised for mismatched chromosome naming"""
```

### Integration Tests
```python
# test_pipeline_atlas_mode.py

def test_full_pipeline_with_test_atlas():
    """End-to-end test with small test BAM and atlas"""
    # Should produce multi-row output for reads with multiple candidates
    
def test_pipeline_v0_single_pas_still_works():
    """Backward compatibility: mode without atlas should work"""
```

---

## 7. References & Related Work

### Relevant Literature
- PolyASite Atlas: [Gruber et al., 2016](https://academic.oup.com/nar/article/44/D1/D86/2502648)
- Alternative polyadenylation: [Mayr, 2016](https://www.nature.com/articles/nrm.2016.51)
- Terminal exon analysis: [Edwalds-Gilbert et al., 1997](https://www.sciencedirect.com/science/article/pii/S0022283697903008)

### Public PAS Atlases Available
- PolyASite (Gruber et al.)
- PolyA-DB
- ENSEMBL PAC (polyadenylation clusters)

---

## 8. Questions for Refinement

1. **Distance threshold:** What value is biologically appropriate? 
   - Default 5 kb seems reasonable for most genes
   - May vary by gene type (long vs short genes)
   - Should be user-configurable

2. **Multiple PAS per gene:** How to handle genes with multiple known PAS?
   - Current design: include all in distance threshold band
   - Alternative: rank by prior evidence, use top N

3. **Unmapped reads:** Reads with no candidate PAS after filtering?
   - Option A: Exclude from output
   - Option B: Include with distance = NA or infinite
   - Recommendation: Option B (report separately for QC)

4. **Performance:** Should we precompute read-to-PAS mapping?
   - Depends on N and M sizes
   - May warrant spatial indexing (e.g., bedtools index)

5. **Output format:** CSV becomes very large (M × N rows)
   - Parquet/HDF5 alternative?
   - Database backend?
   - Multiple output options configurable via CLI?

---

**Last Updated:** 2026-03-04  
**Status:** Ready for review and implementation prioritization
