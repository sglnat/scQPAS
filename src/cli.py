"""
Unified CLI for scQPAS - Single-cell Quantification of PolyAdenylation Sites.

This module provides a single command-line interface that orchestrates the complete
pipeline with in-memory processing (only temporary files required for bedtools).
"""

import click
import pysam
import pandas as pd
import subprocess
import tempfile
import os
from pathlib import Path

from src.reads import extract_reads
from src.annotation import extract_genes, extract_exons, calculate_introns, get_bed_from_df
from src.cigar import get_cigar_bed, filter_by_cigar
from src.distances import calculate_distances


def _run_bedtools_intersect(bed_a_df, bed_b_df, tmpdir, name_a="a", name_b="b", flags=None):
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
        Intersection results
    """
    if flags is None:
        flags = ['-wo']  # default: write overlapping positions from both
    
    bed_a_path = os.path.join(tmpdir, f'bed_{name_a}.bed')
    bed_b_path = os.path.join(tmpdir, f'bed_{name_b}.bed')
    output_path = os.path.join(tmpdir, f'intersection_{name_a}_{name_b}.bed')
    
    # Write temporary BED files
    bed_a_df.to_csv(bed_a_path, sep='\t', header=False, index=False)
    bed_b_df.to_csv(bed_b_path, sep='\t', header=False, index=False)
    
    # Run bedtools intersect with custom flags
    try:
        cmd = ['bedtools', 'intersect'] + flags + ['-a', bed_a_path, '-b', bed_b_path]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        if result.stdout.strip():
            with open(output_path, 'w') as f:
                f.write(result.stdout)
            return pd.read_csv(output_path, sep='\t', header=None)
        else:
            return pd.DataFrame()
    except FileNotFoundError:
        raise click.ClickException(
            "bedtools not found. Please install bedtools: conda install -c bioconda bedtools"
        )
    except subprocess.CalledProcessError as e:
        raise click.ClickException(f"bedtools error: {e.stderr}")





@click.command()
@click.option(
    '--bam',
    type=click.Path(exists=True),
    required=True,
    help='Path to input BAM file (must be indexed)'
)
@click.option(
    '--gtf',
    type=click.Path(exists=True),
    required=True,
    help='Path to GTF annotation file'
)
@click.option(
    '--chr',
    type=str,
    required=True,
    help='Target chromosome (e.g., chr12)'
)
@click.option(
    '--pas',
    type=int,
    required=True,
    help='Polyadenylation site position (genomic coordinate)'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='distances.csv',
    help='Output file path (CSV format). Default: distances.csv'
)
@click.option(
    '--percentage-threshold',
    type=int,
    default=80,
    help='Min %% of A nucleotides in polyA region (0-100). Default: 80'
)
@click.option(
    '--length-threshold',
    type=int,
    default=5,
    help='Min length of polyA tail (bp). Default: 5'
)
@click.option(
    '--use-fc',
    is_flag=True,
    default=False,
    help='Use fixed soft-clipped coordinates from BAM XO/XF tags'
)
def main(bam, gtf, chr, pas, output, percentage_threshold, length_threshold, use_fc):
    """
    Calculate distances from CPA sites to reads.
    
    Processes a BAM file against a GTF annotation in 6 steps:
    1. Extract polyA-containing reads from BAM
    2. Extract annotation (genes, exons, introns) from GTF
    3. Extract CIGAR-derived introns from reads
    4. Intersect reads with annotated features (bedtools)
    5. Filter reads by CIGAR string validation
    6. Calculate distances from cleavage to PAS
    
    All processing in-memory except temporary files for bedtools (auto-cleaned).
    
    Example:
    
        scqpas --bam sample.bam --gtf annotation.gtf --chr chr12 --pas 6538371 --output results.csv
    """
    
    try:
        click.echo("\n" + "=" * 70)
        click.echo("scQPAS: Single-cell Quantification of PolyAdenylation Sites")
        click.echo("=" * 70)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            
            # ========== STEP 1: EXTRACT READS ==========
            click.echo("\n[1/6] Extracting reads from BAM file...")
            with pysam.AlignmentFile(bam, "rb") as sam:
                reads_df = extract_reads(
                    sam,
                    percentage_threshold=percentage_threshold,
                    length_threshold=length_threshold,
                    use_fc=use_fc,
                )
            
            polyA_count = reads_df['is_polyA'].sum()
            click.echo(f"      ✓ Extracted {len(reads_df)} total reads")
            click.echo(f"      ✓ Found {polyA_count} polyA-containing reads")
            
            # ========== STEP 2: EXTRACT ANNOTATION ==========
            click.echo("\n[2/6] Extracting annotation from GTF...")
            
            # All functions now support in-memory operation (output=None)
            genes_df = extract_genes(gtf)
            exons_df = extract_exons(gtf)
            introns_df = calculate_introns(exons_df)
            
            click.echo(f"      ✓ Extracted {len(genes_df)} genes")
            click.echo(f"      ✓ Extracted {len(exons_df)} exons")
            click.echo(f"      ✓ Calculated {len(introns_df)} introns")
            
            # ========== STEP 3: EXTRACT CIGAR-DERIVED INTRONS ==========
            click.echo("\n[3/6] Parsing CIGAR strings...")
            
            # Pass DataFrame directly (no file conversion needed)
            cigar_df = get_cigar_bed(reads_df)
            click.echo(f"      ✓ Extracted {len(cigar_df)} CIGAR-derived introns")
            
            # ========== STEP 4: RUN BEDTOOLS INTERSECTIONS ==========
            click.echo("\n[4/6] Intersecting reads with annotated features...")
            
            # Convert reads and annotation to BED format
            reads_bed_df = get_bed_from_df(reads_df, chr, pas)
            genes_bed_df = genes_df[['chr', 'start', 'end', 'gene_id', 'dummy', 'strand']]
            introns_bed_df = introns_df[['chr', 'start', 'end', 'intron_id', 'length_intron', 'strand']]
            
            # Intersect reads with genes
            # Flags: -wa (write A), -s (same strand), -f 1.0 (A must be 100% covered by B)
            reads_genes_df = _run_bedtools_intersect(
                reads_bed_df, genes_bed_df, tmpdir, 
                name_a="reads", name_b="genes",
                flags=['-wa', '-s', '-f', '1.0']
            )
            click.echo(f"      ✓ Found {len(reads_genes_df)} read-gene intersections")
            
            # Intersect reads with introns
            # Flags: -s (same strand), -F 1.0 (B must be 100% covered by A), -wa -wb (write both)
            reads_introns_df = _run_bedtools_intersect(
                reads_genes_df, introns_bed_df, tmpdir, 
                name_a="reads_genes", name_b="introns",
                flags=['-s', '-F', '1.0', '-wa', '-wb']
            )
            click.echo(f"      ✓ Found {len(reads_introns_df)} read-intron intersections")
            
            # ========== STEP 5: FILTER BY CIGAR ==========
            click.echo("\n[5/6] Filtering reads by CIGAR string validation...")
            
            # Pass DataFrames directly (no file conversion needed)
            valid_reads = filter_by_cigar(
                reads_df,
                reads_genes_df,
                reads_introns_df,
                cigar_df
            )
            
            click.echo(f"      ✓ Retained {len(valid_reads)} valid reads")
            click.echo(f"      ✓ Intronic reads: {valid_reads['intronic'].sum()}")
            
            # ========== STEP 6: CALCULATE DISTANCES ==========
            click.echo("\n[6/6] Calculating distances from CPA sites...")
            
            distances_df = calculate_distances(valid_reads, reads_introns_df)
            
            click.echo(f"      ✓ Calculated distances for {len(distances_df)} reads")
        
        # ========== WRITE OUTPUT ==========
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        distances_df.to_csv(output_path, index=False)
        
        # ========== SUMMARY ==========
        click.echo("\n" + "=" * 70)
        click.echo("✓ PIPELINE COMPLETED SUCCESSFULLY")
        click.echo("=" * 70)
        click.echo(f"Results written to: {output_path.absolute()}")
        click.echo(f"Total reads processed: {len(distances_df)}")
        click.echo(f"Mean distance: {distances_df['distance'].mean():.1f} bp")
        click.echo(f"Median distance: {distances_df['distance'].median():.1f} bp")
        click.echo(f"Std dev: {distances_df['distance'].std():.1f} bp")
        click.echo("=" * 70 + "\n")
        
    except click.ClickException:
        raise
    except Exception as e:
        click.echo(f"\n✗ ERROR: {e}\n", err=True)
        import traceback
        click.echo(traceback.format_exc(), err=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
