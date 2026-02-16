import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from src.distances import df_distances


def main():
    parser = argparse.ArgumentParser(description='Calculate distances from CPA sites to introns')

    # Input arguments
    parser.add_argument('--reads', required=True, help='Path to input reads CSV file')
    parser.add_argument('--valid-reads', required=True, help='Path to valid reads TSV file')
    parser.add_argument('--intronic-reads', required=True, help='Path to intronic reads BED file')

    # Output arguments
    parser.add_argument('--out-distances', required=True, help='Path to output distances CSV file')

    args = parser.parse_args()

    df_distances(args.reads, args.valid_reads, args.intronic_reads, args.out_distances)

if __name__ == '__main__':
    main()


# reads_path = 'GAPDH/TEST_GAPDH_polyA_reads.csv'
# valid_reads_path = 'GAPDH/GAPDH_valid_reads_introns_intersected.tsv'
# intronic_reads_path = 'GAPDH/GAPDH_reads_introns_intersected.bed'
# out_distances_path = 'TEST/TEST_GAPDH_distances9.csv'