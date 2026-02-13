import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from src.annotation import extract_genes, extract_exons, calculate_introns, get_bed_from_df

def main():
    parser = argparse.ArgumentParser(description='Extract annotation from GTF file')
    
    # Input arguments
    parser.add_argument('--gtf', required=True, help='Path to input GTF file')
    parser.add_argument('--reads-csv', required=True, help='Path to reads CSV (for get_bed_from_df)')
    parser.add_argument('--chr', required=True, help='Chromosome to filter (for get_bed_from_df)')
    parser.add_argument('--cpa-site', type=int, required=True, help='CPA site position (for get_bed_from_df)')
    
    # Output arguments
    parser.add_argument('--genes-bed', required=True, help='Path to output genes BED file')
    parser.add_argument('--exons-bed', required=True, help='Path to output exons BED file')
    parser.add_argument('--introns-bed', required=True, help='Path to output introns BED file')
    parser.add_argument('--reads-bed', required=True, help='Path to output reads BED file')
    
    args = parser.parse_args()
    
    # Extract genes
    extract_genes(args.gtf, args.genes_bed)
    
    # Extract exons
    exons = extract_exons(args.gtf, args.exons_bed)
    
    # Calculate introns (uses exons output)
    calculate_introns(exons, args.introns_bed)
    
    # Convert reads CSV to BED
    # this is for distance calculation
    # hardcoded chr and pos will need to be replaced by list (??) of all possible chr and pos for each read!
    get_bed_from_df(args.reads_csv, args.chr, args.cpa_site, args.reads_bed)


if __name__ == '__main__':
    main()

# modify to accept both hardcoded numbers and command line arguments!



# genes = extract_genes('data/v42/gencode.v42.basic.annotation.sorted.gtf', 'GAPDH/gencode.v42.basic.annotation.sorted.gtf.genes.bed')
# exons = extract_exons('data/v42/gencode.v42.basic.annotation.sorted.gtf', 'GAPDH/gencode.v42.basic.annotation.sorted.gtf.exons.bed')
# calculate_introns(exons, 'GAPDH/gencode.v42.basic.annotation.sorted.gtf.introns.bed')

# # this is for distance calculation
# # hardcoded chr and pos will need to be replaced by list (??) of all possible chr and pos for each read!
# get_bed_from_df('GAPDH/TEST_GAPDH_polyA_reads.csv', 'chr12', 6538371, 'GAPDH/TEST_GAPDH_reads_intersect.bed')