import pandas as pd
import numpy as np



def extract_exons(gtf_file_path, exons_output):
    '''
    Extract exons from a GTF file and save them in BED format.
    '''

    gtf_df = pd.read_csv(gtf_file_path, sep='\t', comment='#', header=None)
    # exons_df = gtf_df[gtf_df[2].isin(['transcript', 'exon'])]
    exons = gtf_df[gtf_df[2] == 'exon']
    exons = exons[[0, 3, 4, 8, 6]].copy() # 0: chromosome, 3: start, 4: end, 6: strand, 8: attributes
    exons.columns = ['chr', 'start', 'end', 'attributes', 'strand']
    exons['start'] = exons['start'] - 1  # Convert start to 0-based
    exons['transcript_id'] = exons['attributes'].str.extract(r'transcript_id\s+"([^"]+)"')
    exons['gene_id'] = exons['attributes'].str.extract(r'gene_id\s+"([^"]+)"')
    exons.drop(columns='attributes', inplace=True)
    exons = exons[['chr', 'start', 'end', 'transcript_id', 'strand', 'gene_id']]
    #print(exons)
    exons.to_csv(exons_output, sep='\t', header=False, index=False, columns=['chr', 'start', 'end', 'transcript_id', 'strand', 'gene_id'])

    return exons



def calculate_introns(exons, introns_output):
    '''
    Returns dataframe with transcript-specific introns, defined as the regions between exons of the same transcript
    '''

    grouped_exons = exons.groupby('transcript_id')

    introns_list = []

    for transcript_id, group in grouped_exons:
        group = group.sort_values(by='start')
        exons_starts = group['start'].to_list()
        exons_ends = group['end'].to_list()

        # Calculate intron positions
        for i in range(len(exons_starts) - 1):
            intron_start = exons_ends[i] # points to the first bp of the intron
            intron_end = exons_starts[i + 1] # points one bp past the last bp of the intron, as per .bed specifications
            
            #if intron_start <= intron_end:
            intron_length = intron_end - intron_start
            introns_list.append({
                'chr': group['chr'].iloc[0],
                'start': intron_start,
                'end': intron_end,
                'intron_id': f"{transcript_id}_{i}",  # Unique intron identifier
                'length_intron': intron_length,
                'strand': group['strand'].iloc[0],
                'gene_id': group['gene_id'].iloc[0]
            })

    # Create a DataFrame for introns and save as BED file
    introns_df = pd.DataFrame(introns_list) #columns: chr, start, end, intron_id, length_intron, strand
    #print(introns_df)

    # Save to BED format (requires specific format: no header and only certain columns)
    introns_df.to_csv(introns_output, sep='\t', header=False, index=False)



def get_bed_from_df(df_input, chr, cpa_site, bed_output):
    '''
    Save a DataFrame in BED format.
    Expects the DataFrame to have columns: 'chr', 'start', 'end', 'name', 'strand', 'gene_id'
    '''

    #df = pd.read_csv(df_input, sep=',', names=['UMI', 'CB', 'chr', 'start', 'end', 'strand', 'CIGAR', 'is_polyA', 'len_pA', 'cpa_site', 'is_polyA_RS', 'distance_to_cpa'])
    df = pd.read_csv(df_input, sep=',')

    print(df.head())

    bed = []
    
    for idx, row in df.iterrows():
        if row['strand'] == '+' and row['chr'] == chr:
            start = row['start']
            end = cpa_site + 1  # BED end is exclusive
        elif row['strand'] == '-' and row['chr'] == chr:
            start = cpa_site  # BED start is inclusive
            end = row['end']
        else:
            continue  # Skip rows with an undefined strand

        if start < end:
            bed.append([row['chr'], start, end, row['read_id'], row['strand']])

    bed_df = pd.DataFrame(bed, columns=['chr', 'start', 'end', 'read_id', 'strand'])
    bed_df.insert(4, 'dummy', 1000)
    print(bed_df)

    bed_df.to_csv(bed_output, sep='\t', header=False, index=False)



def extract_genes(gtf_file_path, genes_output):
    '''
    Extract genes from a GTF file and save them in BED format.
    '''

    gtf_df = pd.read_csv(gtf_file_path, sep='\t', comment='#', header=None)
    # exons_df = gtf_df[gtf_df[2].isin(['transcript', 'exon'])]
    genes = gtf_df[gtf_df[2] == 'gene']
    genes = genes[[0, 3, 4, 8, 6]].copy() # 0: chromosome, 3: start, 4: end, 6: strand, 8: attributes
    genes.columns = ['chr', 'start', 'end', 'attributes', 'strand']
    genes['start'] = genes['start'] - 1  # Convert start to 0-based
    genes['gene_id'] = genes['attributes'].str.extract(r'gene_id\s+"([^"]+)"')
    genes.drop(columns='attributes', inplace=True)
    genes = genes[['chr', 'start', 'end', 'gene_id', 'strand']]
    genes.insert(4, 'dummy', 1000)
    print(genes)

    genes.to_csv(genes_output, sep='\t', header=False, index=False, columns=['chr', 'start', 'end', 'gene_id', 'dummy', 'strand'])

    #return genes