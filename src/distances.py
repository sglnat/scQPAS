import pandas as pd
import numpy as np



def df_distances(reads_path, valid_reads_path, intronic_reads_path, out_distances_path):
    df = pd.read_csv(reads_path, index_col=0)
    #print(df)
    valid_reads = pd.read_csv(valid_reads_path, sep='\t')
    #print(valid_reads)
    intronic_reads = pd.read_csv(intronic_reads_path, sep='\t', header=None, names=['chr_read', 'start_read', 'end_read', 'read_id', 'score_read', 'strand_read', 'chr_intron', 'start_intron', 'end_intron', 'intron_id', 'length_intron', 'strand_intron', 'gene_id'])
    #print(intronic_reads)
    # Determine the best CPA site
    df_cpa_sites = get_cpa_sites(df)
    #print(df_cpa_sites)
    #best_cpa_site = best_cpa(df_cpa_sites)
    #print('Best CPA site: ' + str(best_cpa_site))
    df_distance = calculate_distances(valid_reads, intronic_reads)
    df_distance.to_csv(out_distances_path, index=False)

####### Modify so that it takes multiple CPA sites #######



def get_cpa_sites(df):    
    """
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing read information including cpa_site column.

    Returns
    -------            
    df_cpa_sites : pd.DataFrame
        Updated DataFrame with new cpa sites.

    """       

    # Initialize a DataFrame to hold the CPA site counts
    df_cpa_sites = pd.DataFrame(columns=['cpa_site', 'supporting_read_count'])
    
    # Count occurrences of each CPA site in the df
    counts = df['cpa_site'].value_counts()
    
    # Create a DataFrame from the counts Series
    df_cpa_sites = counts.reset_index()
    df_cpa_sites.columns = ['cpa_site', 'supporting_read_count']

    return df_cpa_sites



def best_cpa(df_cpa_sites):    
    best_read_count = 0
    best_cpa_site = None

    for _, row in df_cpa_sites.iterrows():
        cpa_site = row['cpa_site']
        read_count = row[1]
    
        if read_count > best_read_count:
            best_read_count = read_count
            best_cpa_site = cpa_site

    return best_cpa_site



def calculate_distances(valid_reads, intronic_reads):

    intronic_reads['transcript_id'] = intronic_reads['intron_id'].str.split('_').str[0]
    #print(intronic_reads)

    # Create a new DataFrame for storing distances
    distance_df = pd.DataFrame()
    distance_df[['read_id', 'transcript_id', 'intronic', 'CIGAR_N']] = valid_reads[['read_id', 'transcript_id', 'intronic', 'CIGAR_N']]
    #print(distance_df)

    # Merging to include total intron lengths in the distance_df
    distance_df = distance_df.merge(intronic_reads.groupby(['read_id', 'transcript_id'])['length_intron'].sum().reset_index(), on=['read_id', 'transcript_id'], how='left')
    distance_df.rename(columns={'length_intron': 'total_intron_len'}, inplace=True)

    # Calculate distance based on strand information
    distance_df['distance'] = None

    # Calculate distance for strand '+'
    mask_plus = valid_reads['strand'] == '+'
    distance_df.loc[mask_plus, 'distance'] = valid_reads.loc[mask_plus, 'end'].values - valid_reads.loc[mask_plus, 'start'].values

    # Calculate distance for strand '-'
    mask_minus = valid_reads['strand'] == '-'
    distance_df.loc[mask_minus, 'distance'] = valid_reads.loc[mask_minus, 'start'].values - valid_reads.loc[mask_minus, 'end'].values

    # Subtract total_intron_len from distance for intronic reads
    intronic_mask = valid_reads['intronic'] == True
    distance_df.loc[intronic_mask, 'distance'] -= distance_df.loc[intronic_mask, 'total_intron_len']

    #print(distance_df)

    return distance_df