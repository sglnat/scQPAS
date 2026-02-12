import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from src.gtf_extend import get_adj_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Adjust terminal exons in GTF file.')
    parser.add_argument('input_file', type=str, nargs='?', default='../data/v42/gencode.v42.basic.annotation.sorted.gtf', help='Full path to the input GTF file')
    parser.add_argument('output_file', type=str, nargs='?', default='../TEST/gencode.v42.basic.annotation.sorted_TEST_pipeline.gtf', help='Full path for the modified GTF file')

    args = parser.parse_args()

    get_adj_df(args.input_file, args.output_file)