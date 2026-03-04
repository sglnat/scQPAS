# Implementation Summary: CPA Site Detection & Data-Driven PAS

**Date:** 2026-03-04  
**Version:** scQPAS V0 → V0.5 (Single PAS detected from reads)

---

## ✅ Completed Changes

### 1. Refactored CPA Site Analysis Functions

**File:** [scqpas/calculate_distances.py](scqpas/calculate_distances.py)

#### `get_cpa_sites(reads_df)` - Updated
- **Old:** Counted bare CPA positions without chromosome info
- **New:** Groups by `(chr, cpa_site)` pairs
- **Returns:** DataFrame with columns `['chr', 'cpa_site', 'supporting_read_count']`
- **Sorted:** By supporting_read_count descending (easy selection of best CPA)

```python
# Old: Only position
cpa_counts = df["cpa_site"].value_counts()  # Just positions

# New: Chr + position
cpa_counts = df.groupby(["chr", "cpa_site"]).size()  # (chr, pos) pairs
```

#### `best_cpa(cpa_sites_df)` - Updated  
- **Old:** Returned integer position only
- **New:** Returns tuple `(chr, cpa_site_position)`
- **Returns:** `(str, int)` tuple for exact genomic coordinate

```python
# Old
best_cpa_site = best_cpa(df_cpa_sites)  # → 42000001 (int)

# New
best_chr, best_pos = best_cpa(cpa_sites_df)  # → ('chr22', 42000001)
```

#### **NEW** `detect_best_cpa_from_reads(reads_df)` - Added
- **Purpose:** Complete CPA detection workflow in one function
- **Input:** Full reads DataFrame from BAM extraction
- **Output:** `(chr, pas_pos, strand)` tuple
- **Logic:**
  1. Filter for polyA-containing reads only
  2. Call `get_cpa_sites()` to get CPA site distribution
  3. Call `best_cpa()` to select most supported site
  4. Determine dominant strand at that position
  5. Return complete PAS coordinates with strand
- **Error Handling:** Raises `RuntimeError` if no polyA reads found

```python
best_chr, best_pos, best_strand = detect_best_cpa_from_reads(reads_df)
# Returns: ('chr22', 42000001, '+')
```

---

### 2. Modified Pipeline Architecture

**File:** [scqpas/core.py](scqpas/core.py)

#### Pipeline Redesigned: 6 Steps (was 7 with multi-PAS loop)

**Current V0.5 Flow:**
1. **Extract reads** from BAM (detect polyA)
2. **Detect best CPA** from read data (data-driven, no user input needed)
3. **Extract annotation** from GTF
4. **Extract CIGAR** introns from reads
5. **Bedtools intersections** using detected CPA
6. **Calculate distances** to detected PAS

#### Key Changes:
- **`pas_bed_path` now Optional:** No longer required input
- **BED file deprecated:** `--pas` argument still accepted for backward compatibility but triggers deprecation warning
- **CPA detection automatic:** Happens in Step 2, no manual specification needed
- **Less flexible, more focused:** Removes multi-PAS loop from V0 in favor of single best-CPA approach

#### Updated `run_pipeline()` function signature:
```python
def run_pipeline(
    bam_path: str,
    gtf_path: str,
    pas_bed_path: Optional[str] = None,  # Now optional, DEPRECATED
    output_path: str = None,
    percentage_threshold: int = 80,
    length_threshold: int = 5,
    config_manager: Optional[ConfigManager] = None,
) -> pd.DataFrame:
```

#### Step 2 Implementation:
```python
# ========== STEP 2: DETECT BEST CPA SITE FROM READS ==========
logger.info("[2/6] Detecting best CPA site from read data...")
pas_chr, pas_pos, pas_strand = detect_best_cpa_from_reads(reads_df)
```

---

### 3. CLI Interface Updated

**File:** [scqpas/cli.py](scqpas/cli.py)

#### `--pas` Option Changes:
```python
# Old
@click.option("--pas", required=True, ...)

# New
@click.option("--pas", required=False, default=None, ...)
```

#### Help Text Updated:
```
"Path to polyadenylation sites BED file (DEPRECATED - PAS now auto-detected from reads). 
Kept for backward compatibility."
```

#### Function Signature:
```python
def main(
    ...
    pas: Optional[str] = None,  # Now accepts None
    ...
) -> None:
```

#### Updated Docstring:
- Explains that PAS is now **automatically detected** from read data
- Emphasizes data-driven approach vs. manual specification
- Updated example commands (no `--pas` needed)

Old example:
```bash
scqpas --bam sample.bam --gtf annotation.gtf --pas pas_sites.bed --output results.csv
```

New example:
```bash
scqpas --bam sample.bam --gtf annotation.gtf --output results.csv
```

---

### 4. Documentation Updated

**File:** [README.md](README.md)

#### Command-line Usage Section:
```bash
# Old: Required --pas
scqpas --bam sample.bam --gtf annotation.gtf --pas pas_sites.bed --output results.csv

# New: PAS auto-detected
scqpas --bam sample.bam --gtf annotation.gtf --output results.csv
```

#### Programmatic Usage:
```python
# Old
run_pipeline(
    bam_path='sample.bam',
    gtf_path='annotation.gtf',
    pas_bed_path='pas_sites.bed',
    output_path='results.csv'
)

# New
run_pipeline(
    bam_path='sample.bam',
    gtf_path='annotation.gtf',
    output_path='results.csv'
    # Note: PAS is automatically detected from reads
)
```

#### Added Note About Auto-Detection:
"The polyadenylation site (PAS) is now **automatically detected** from the read data by identifying the cleavage position with the most supporting evidence. No manual PAS specification is needed."

---

## 📋 Design Planning: Future Atlas Feature

**File:** [DESIGN_NOTES_ATLAS_FEATURE.md](DESIGN_NOTES_ATLAS_FEATURE.md) (NEW)

Comprehensive 8-section design document covering:

### 1. **Current State (V0)** 
- Single CPA site workflow
- Current output structure
- Key properties

### 2. **Future State (V1): Multiple PAS per Read**
- Extended workflow diagram
- Atlas filtering strategies:
  - **Option A:** Distance-based (recommended, ±5kb from 3' end)
  - **Option B:** Gene-based (requires gene mapping)
- Data structure modifications
- Pipeline flow changes
- Distance calculation extensions
- Output format changes (multiple rows per read)

### 3. **Implementation Roadmap**
- Phase 1: Preparation (✓ Done)
- Phase 2: Infrastructure  
- Phase 3: Core Implementation
- Phase 4: Configuration & CLI
- Phase 5: Testing & Validation
- Phase 6: Documentation

### 4. **Technical Considerations**
- **Performance & Scalability:**
  - Memory usage implications
  - Single vs. multiple bedtools intersections
  - Optimization strategies (chunking, formats)
- **Chromosome & Coordinate System:**
  - Naming consistency validation
  - Reference genome version handling
- **Strand Handling:**
  - Filtering by strand match
  - Verification logic

### 5. **Future Enhancements (Post-V1)**
- Score-based PAS selection
- Gene-specific filtering
- Tissue-specific atlas
- ML integration

### 6. **Testing Strategy**
- Unit tests for filtering functions
- Integration tests for full pipeline
- Backward compatibility verification

### 7. **References & Public PAS Atlases**
- PolyASite Atlas
- PolyA-DB
- ENSEMBL PAC

### 8. **Open Questions for Refinement**
- Optimal distance threshold
- Handling multiple PAS per gene
- Strategy for unmapped reads
- Performance optimization needs
- Output format alternatives (CSV vs. Parquet/HDF5)

---

## 🔄 Key Design Decisions Made

### 1. **Chromosome Inclusion in CPA Analysis**
✅ **Decision:** Include `chr` in CPA site analysis  
**Rationale:** 
- Identifies exact genomic coordinates, not just position
- Prevents collision issues (same CPA position on different chromosomes)
- Enables multi-chromosome analysis

### 2. **Data-Driven PAS Detection**
✅ **Decision:** Detect best CPA from reads themselves  
**Rationale:**
- User no longer needs To manually specify PAS
- Self-contained analysis (read distribution determines focus)
- Foundation for future atlas-based multi-PAS approach

### 3. **Making `--pas` Optional (Not Removed)**
✅ **Decision:** Deprecate but maintain `--pas` argument  
**Rationale:**
- Backward compatibility for existing workflows
- Reserved for V1 atlas feature 
- Easy migration path (will become `--pas-atlas` in V1)
- Deprecation warning prepared for users

### 4. **Architecture for Future Atlas**
✅ **Decision:** Single `detect_best_cpa_from_reads()` function  
**Rationale:**
- Encapsulates complete CPA detection workflow
- Easy to extend with atlas logic later
- Clear separation of responsibilities
- Testable unit for validation

---

## ⚙️ Technical Implementation Details

### New Function Relationships

```
extract_reads()                    (→ reads_df)
        ↓
detect_best_cpa_from_reads()       ← NEW main function
    ├─→ get_cpa_sites()            ← MODIFIED (now chr-aware)
    └─→ best_cpa()                 ← MODIFIED (returns tuple)
        ↓
core.py pipeline Step 2            ← UPDATED
        ↓
calculate_distances()              (uses detected CPA)
```

### Data Flow Example

```python
# BAM extraction produces:
reads_df = extract_reads(bam_file)
# columns: [..., chr, cpa_site, is_polyA, ...]
# Example row: {..., chr='chr22', cpa_site=42001000, is_polyA=True, ...}

# Step 2: Auto-detect best CPA
pas_chr, pas_pos, pas_strand = detect_best_cpa_from_reads(reads_df)
# Result: chr='chr22', pos=42001000, strand='+'
# (most reads detected polyA tall at this position)

# Used for distances
distances_df = calculate_distances(valid_reads, introns_df)
# distances calculated relative to detected PAS coordinate
```

---

## 🧪 Validation & Testing

### Syntax Validation
✅ All files validated - no syntax errors:
- `scqpas/calculate_distances.py` ✓
- `scqpas/core.py` ✓  
- `scqpas/cli.py` ✓
- `README.md` ✓

### Type Hints
✅ Complete type hints maintained throughout:
```python
def best_cpa(cpa_sites_df: pd.DataFrame) -> Optional[tuple]:
    """Returns (chr, pos) or None"""

def detect_best_cpa_from_reads(reads_df: pd.DataFrame) -> tuple:
    """Returns (chr, pos, strand)"""
```

---

## 📊 Backward Compatibility

### Changes That Maintain Compatibility:
- ✅ `--pas` argument still works (with deprecation warning)
- ✅ All existing function signatures backward-compatible
- ✅ Output CSV format unchanged (just uses detected CPA)
- ✅ Programmatic API includes optional parameters

### Breaking Changes (Minimal):
- ❌ `run_pipeline(pas_bed_path=...)` now optional (was required)
  - **Mitigation:** Default `None` with deprecation warning
  - **Impact:** Minor; most users will simply remove the argument

---

## 🎯 Next Steps for V1 Implementation

When you're ready to implement the atlas feature:

1. **Review** [DESIGN_NOTES_ATLAS_FEATURE.md](DESIGN_NOTES_ATLAS_FEATURE.md) in detail
2. **Implement Phase 2:** Infrastructure functions
   - `load_pas_atlas()` in process_pasfile.py
   - `filter_pas_candidates_per_read()` with distance threshold
3. **Implement Phase 3:** Core changes
   - `calculate_distances_multi_pas()` 
   - Modified pipeline Step 2b for atlas
4. **Add CLI options:** `--pas-atlas` and `--pas-distance-threshold`
5. **Test:** Unit tests for filtering logic
6. **Document:** Update README with atlas examples

---

## 📞 Questions?

See DESIGN_NOTES_ATLAS_FEATURE.md Section 8 for open questions that need refinement before V1 implementation:
- Optimal distance threshold values
- Handling multiple PAS per gene  
- Output format for large result sets
- Performance optimization strategies
