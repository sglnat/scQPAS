"""
Core pipeline logic for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module contains the main pipeline orchestration logic, independent of CLI interface.
All processing is in-memory except temporary files required for bedtools (auto-cleaned).
"""

import logging
import tempfile
from pathlib import Path
import pysam

from .extract_reads_BAM import extract_reads
from .extract_annotation_GTF import extract_genes, extract_exons, calculate_introns, get_bed_from_df
from .extract_cigar_BAM import get_cigar_bed, filter_by_cigar
from .calculate_distances import calculate_distances
from .bedtools_intersections import run_bedtools_intersect

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path,
    gtf_path,
    chr,
    pas,
    output_path,
    percentage_threshold=80,
    length_threshold=5,
    use_fc=False,
    config_manager=None
):
    """
    Execute the complete scQPAS pipeline.
    
    Processes a BAM file against a GTF annotation in 6 steps:
    1. Extract polyA-containing reads from BAM
    2. Extract annotation (genes, exons, introns) from GTF
    3. Extract CIGAR-derived introns from reads
    4. Intersect reads with annotated features (bedtools)
    5. Filter reads by CIGAR string validation
    6. Calculate distances from cleavage to PAS
    
    Args:
        bam_path: Path to input BAM file (must be indexed)
        gtf_path: Path to GTF annotation file
        chr: Target chromosome (e.g., 'chr12')
        pas: Polyadenylation site position (genomic coordinate)
        output_path: Output file path for results CSV
        percentage_threshold: Min % of A nucleotides in polyA region (0-100). Default: 80
        length_threshold: Min length of polyA tail (bp). Default: 5
        use_fc: Use fixed soft-clipped coordinates from BAM XO/XF tags. Default: False
        config_manager: ConfigManager instance for accessing pipeline configuration. Default: None
    
    Returns:
        distances_df: DataFrame with calculated distances
    
    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If pipeline processing fails
    """
    
    logger.info("=" * 70)
    logger.info("scQPAS: Single-cell Quantification of PolyAdenylation Sites")
    logger.info("=" * 70)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        
        # ========== STEP 1: EXTRACT READS ==========
        logger.info("[1/6] Extracting reads from BAM file...")
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            reads_df = extract_reads(
                sam,
                percentage_threshold=percentage_threshold,
                length_threshold=length_threshold,
                use_fc=use_fc,
                config_manager=config_manager
            )
        
        polyA_count = reads_df['is_polyA'].sum()
        logger.info(f"      ✓ Extracted {len(reads_df)} total reads")
        logger.info(f"      ✓ Found {polyA_count} polyA-containing reads")
        
        # ========== STEP 2: EXTRACT ANNOTATION ==========
        logger.info("[2/6] Extracting annotation from GTF...")
        
        # All functions now support in-memory operation (output=None)
        genes_df = extract_genes(gtf_path, config_manager=config_manager)
        exons_df = extract_exons(gtf_path, config_manager=config_manager)
        introns_df = calculate_introns(exons_df, config_manager=config_manager)
        
        logger.info(f"      ✓ Extracted {len(genes_df)} genes")
        logger.info(f"      ✓ Extracted {len(exons_df)} exons")
        logger.info(f"      ✓ Calculated {len(introns_df)} introns")
        
        # ========== STEP 3: EXTRACT CIGAR-DERIVED INTRONS ==========
        logger.info("[3/6] Parsing CIGAR strings...")
        
        # Pass DataFrame directly (no file conversion needed)
        cigar_df = get_cigar_bed(reads_df)
        logger.info(f"      ✓ Extracted {len(cigar_df)} CIGAR-derived introns")
        
        # ========== STEP 4: RUN BEDTOOLS INTERSECTIONS ==========
        logger.info("[4/6] Intersecting reads with annotated features...")
        
        # Convert reads and annotation to BED format
        reads_bed_df = get_bed_from_df(reads_df, chr, pas, config_manager=config_manager)
        genes_bed_df = genes_df[['chr', 'start', 'end', 'gene_id', 'dummy', 'strand']]
        introns_bed_df = introns_df[['chr', 'start', 'end', 'intron_id', 'length_intron', 'strand']]
        
        # Intersect reads with genes
        # Flags: -wa (write A), -s (same strand), -f 1.0 (A must be 100% covered by B)
        reads_genes_df = run_bedtools_intersect(
            reads_bed_df, genes_bed_df, tmpdir, 
            name_a="reads", name_b="genes",
            flags=['-wa', '-s', '-f', '1.0']
        )
        
        # Assign column names: -wa outputs columns from reads_bed (6 cols)
        if not reads_genes_df.empty:
            reads_genes_df.columns = ['chr', 'start', 'end', 'read_id', 'dummy', 'strand']
        
        logger.info(f"      ✓ Found {len(reads_genes_df)} read-gene intersections")
        
        # Intersect reads with introns
        # Flags: -s (same strand), -F 1.0 (B must be 100% covered by A), -wa -wb (write both)
        reads_introns_df = run_bedtools_intersect(
            reads_genes_df, introns_bed_df, tmpdir, 
            name_a="reads_genes", name_b="introns",
            flags=['-s', '-F', '1.0', '-wa', '-wb']
        )
        
        # Assign column names: -wa outputs A columns, -wb outputs B columns
        if not reads_introns_df.empty:
            reads_introns_df.columns = [
                'chr_read', 'start_read', 'end_read', 'read_id', 'dummy_read', 'strand_read',
                'chr_intron', 'start_intron', 'end_intron', 'intron_id', 'length_intron', 'strand_intron'
            ]
        
        logger.info(f"      ✓ Found {len(reads_introns_df)} read-intron intersections")
        
        # ========== STEP 5: FILTER BY CIGAR ==========
        logger.info("[5/6] Filtering reads by CIGAR string validation...")
        
        # Pass DataFrames directly (no file conversion needed)
        valid_reads = filter_by_cigar(
            reads_df,
            reads_genes_df,
            reads_introns_df,
            cigar_df
        )
        
        logger.info(f"      ✓ Retained {len(valid_reads)} valid reads")
        logger.info(f"      ✓ Intronic reads: {valid_reads['intronic'].sum()}")
        
        # ========== STEP 6: CALCULATE DISTANCES ==========
        logger.info("[6/6] Calculating distances from CPA sites...")
        
        distances_df = calculate_distances(valid_reads, reads_introns_df)
        
        logger.info(f"      ✓ Calculated distances for {len(distances_df)} reads")
    
    # ========== WRITE OUTPUT ==========
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    distances_df.to_csv(output_path_obj, index=False)
    
    logger.info("=" * 70)
    logger.info("✓ PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 70)
    logger.info(f"Results written to: {output_path_obj.absolute()}")
    logger.info(f"Total reads processed: {len(distances_df)}")
    logger.info("=" * 70)
    
    return distances_df
