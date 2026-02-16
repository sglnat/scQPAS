import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from src.cigar import get_cigar_bed, filter_by_cigar

def main():
    parser = argparse.ArgumentParser(description='Filter reads based on whether their CIGAR string matches the introns from at least one transcript from the annotation')
    
    # Input arguments
    parser.add_argument('--reads-csv', required=True, help='Path to reads CSV (for get_bed_from_df)')
    parser.add_argument('--reads-in-genes', required=True, help='Path to reads within single gene context BED file')
    parser.add_argument('--reads-in-introns', required=True, help='Path to intronic reads BED file')

    # Output arguments
    parser.add_argument('--out-cigar', required=True, help='Path to output CIGAR-based introns')
    parser.add_argument('--out-reads', required=True, help='Path to output CIGAR-validated reads')
    
    args = parser.parse_args()
    
    # Extract CIGAR-based introns
    cf = get_cigar_bed(args.reads_csv, args.out_cigar)
    
    # Filter by CIGAR-based introns
    filter_by_cigar(args.reads_csv, args.reads_in_genes, args.reads_in_introns, cf, args.out_reads)

if __name__ == "__main__":
    main()


# cf = get_cigar_bed('GAPDH/TEST_GAPDH_polyA_reads.csv', 'GAPDH/TEST_GAPDH_reads_cigar.csv')
# filter_by_cigar('GAPDH/TEST_GAPDH_polyA_reads.csv', 'GAPDH/GAPDH_reads_genes_intersected.bed', 'GAPDH/GAPDH_reads_introns_intersected.bed', cf, 'GAPDH/GAPDH_valid_reads_introns_intersected.tsv')
