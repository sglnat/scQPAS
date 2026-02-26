import pandas as pd
import csv

def adjust_terminal_exons(df):
    for gene_id, gene_group in df.groupby('gene_id'):
        for transcript_id, transcript_group in gene_group.groupby('transcript_id'):
            if 'exon' in transcript_group['feature'].values:
                gene_idx = gene_group.index[gene_group['feature'] == 'gene'][0]
                transcript_idx = transcript_group.index[transcript_group['feature'] == 'transcript'][0]
                last_exon_idx = transcript_group[transcript_group['feature'] == 'exon'].index[-1]
                strand = transcript_group.at[last_exon_idx, 'strand']

                if strand == '+':
                    df.at[last_exon_idx, 'end'] += 1000
                    df.at[transcript_idx, 'end'] += 1000
                    df.at[gene_idx, 'end'] += 1000
                    
                elif strand == '-':
                    df.at[last_exon_idx, 'start'] -= 1000
                    df.at[transcript_idx, 'start'] -= 1000
                    df.at[gene_idx, 'start'] -= 1000
                    
    return df

def get_adj_df(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', comment='#', header=None)
    df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    #Extract transcript_id and gene_id for grouping
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)"')
    df['gene_id'] = df['attribute'].str.extract(r'gene_id "([^"]+)"')

    adjusted_df = adjust_terminal_exons(df)

    #Drop temporary columns
    adjusted_df = adjusted_df.drop(columns=['transcript_id', 'gene_id'])

    adjusted_df.to_csv(output_file, sep='\t', header=False, index=None, quoting=csv.QUOTE_NONE)
    #print(f"Adjusted GTF file saved as: {output_file}")