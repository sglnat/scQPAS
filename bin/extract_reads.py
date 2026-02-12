import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')


import pysam
from src.reads import extract_reads

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--percentage-threshold", type=int, default=80)
    parser.add_argument("--length-threshold", type=int, default=5)
    parser.add_argument("--use-fc", action="store_true")
    args = parser.parse_args()

    sam = pysam.AlignmentFile(args.bam, "rb")
    df = extract_reads(
        sam,
        percentage_threshold=args.percentage_threshold,
        length_threshold=args.length_threshold,
        use_fc=args.use_fc,
    )
    df.to_csv(args.output, index=False)

    sam.close()

if __name__ == "__main__":
    main()


# original code for running:
# bam_file = pysam.AlignmentFile("GAPDH/full_GAPDH.bam", "rb")
# polyA_df = extract_reads(bam_file, percentage_threshold=80, length_threshold=5, use_fc=False)
# polyA_df.to_csv("GAPDH/TEST_GAPDH_polyA_reads.csv", index=True)