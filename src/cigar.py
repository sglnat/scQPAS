import pandas as pd
import numpy as np
import re



def get_cigar_bed(df_input, bed_output):
    '''
    Save a DataFrame with CIGAR strings in BED format.
    Expects the DataFrame to have columns: 'chr', 'start', 'end', 'name', 'strand', 'CIGAR'
    '''

    #unmatched_df = pd.DataFrame(columns=['chr', 'start', 'end', 'read_id', 'strand', 'CIGAR'])
    unmatched_list = []

    # Split the CIGAR string into segments
    df = pd.read_csv(df_input, sep=',')

    print(df)
    
    for _, row in df.iterrows():
        cig_segments = re.findall(r'(\d+)([MIDNSHPX=])', row['CIGAR'])

        # Loop through segments to find unmatched N segments
        # Since CIGAR is always on + strand, this calculation works for both + and - strand
        current_start = row['start']
        for length, operation in cig_segments:
                length = int(length) 
                if operation == 'N':
                    unmatched_row = {
                        'chr': row['chr'],
                        'start': current_start,
                        'end': current_start + length,
                        'read_id': row['read_id'],
                        'strand': row['strand'],
                        #'CIGAR_seq': f'{length}N'
                        'CIGAR': row['CIGAR']
                    }
                    unmatched_list.append(unmatched_row)
                elif operation in ['M', 'D', '=', 'X']:
                    current_start += length

    unmatched_df = pd.DataFrame(unmatched_list)

    print(unmatched_df)
    unmatched_df.to_csv(bed_output, sep='\t', header=False, index=False)

    return unmatched_df


def filter_by_cigar(reads, filtered_reads, intronic_reads, cigar_df, output_csv):
    '''
    - Keeps reads where all introns obtained by CIGAR string match the introns from at least one transcript from the annotation. 
    If they do not match any transcript, they are from an unannotated transcript isoform, which cannot be accounted for in the distance calculation, so they are excluded.
    - Keeps reads without any introns in CIGAR string (no N operation) where the distance to the PAS covers at least one intron.
    - Keeps reads close to the cleavage site where the distance to the PAS does not cover any intron.
    - By using filtered_reads, reads outside of a single gene context are excluded.

    These four cases should cover all reads.
    '''

    reads = pd.read_csv(reads, sep=',')
    filtered_reads = pd.read_csv(filtered_reads, sep='\t', header=None, names=['chr', 'start', 'end', 'read_id', 'dummy', 'strand'])
    cigar_df.columns = ['chr_cf', 'start_cf', 'end_cf', 'read_id', 'strand_cf', 'CIGAR']
    intronic_reads = pd.read_csv(intronic_reads, sep='\t', header=None, names=['chr_read', 'start_read', 'end_read', 'read_id', 'score_read', 'strand_read', 'chr_intron', 'start_intron', 'end_intron', 'intron_id', 'length_intron', 'strand_intron', 'gene_id'])
    intronic_reads['transcript_id'] = intronic_reads['intron_id'].str.split('_').str[0]

    # Merge DataFrames on read_id
    merged_df = pd.merge(cigar_df, intronic_reads, on='read_id', suffixes=('_cf', '_intron'), how='left')
    # Create a boolean series for matches
    merged_df['match'] = (merged_df['start_cf'] == merged_df['start_intron']) & (merged_df['end_cf'] == merged_df['end_intron'])
    #merged_df.to_csv('TEST/merged_df.tsv', sep='\t', index=False)

    # Group by read_id and transcript_id to check if all introns matched
    intronic_N = merged_df.groupby(['read_id', 'transcript_id']).agg(all_match=('match', 'all')).reset_index().query('all_match').drop(columns='all_match')

    # Keep reads from df that do not have an N in their CIGAR string (i.e., no introns)
    cigar_noN = reads[reads['CIGAR'].str.contains('N') == False]['read_id']
    intronic_noN = intronic_reads[intronic_reads['read_id'].isin(cigar_noN)][['read_id', 'transcript_id']].groupby(['read_id', 'transcript_id'], as_index=False).first()

    # Mark them as intronic
    intronic_N['intronic'] = True
    intronic_N['CIGAR_N'] = True
    intronic_noN['intronic'] = True
    intronic_noN['CIGAR_N'] = False

    # Handle reads close to the cleavage site -> excluding those already investigated as intronic
    excluded_reads = pd.concat([intronic_reads[['read_id']], cigar_df[['read_id']]])
    close_reads_list = filtered_reads['read_id']
    close_reads = close_reads_list[~close_reads_list.isin(excluded_reads['read_id'])]

    close_reads = pd.DataFrame({
        'read_id': close_reads,
        'transcript_id': np.nan,  
        'intronic': False,
        'CIGAR_N': False          
    }) 

    # Combine intronic reads with valid CIGAR N, reads without CIGAR N, and non-intronic reads
    final_reads = pd.concat([intronic_N, intronic_noN, close_reads], ignore_index=True)
    # Merge final_reads with filtered_reads to keep chr, start, end
    final_reads = pd.merge(final_reads, filtered_reads[['read_id', 'chr', 'start', 'end', 'strand']], on='read_id', how='left')
    final_reads.to_csv(output_csv, sep='\t', index=False)

    return final_reads