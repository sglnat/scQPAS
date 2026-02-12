import pandas as pd
import numpy as np

def extract_reads(sam, percentage_threshold, length_threshold, use_fc):
    reads = []
    read_ids = []
    
    for read in sam.fetch():
        read_id = read.query_name
        umis = read.get_tag('UB') if read.has_tag('UB') else None 
        cb = read.get_tag('CB') if read.has_tag('CB') else None 
        chr_name = read.reference_name
        start = read.reference_start # 0-based leftmost coordinate
        end = read.reference_end # points one past the last aligned base
        rev = read.is_reverse
        strand = '-' if rev == True else '+' 
        cigar = read.cigarstring

        # CIGAR tuples for softclip
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]

        # Check if the read is polyA
        is_polyA, len_pA = check_polyA(read, left_end, right_end, rev, percentage_threshold, length_threshold, use_fc)

        # Update the CPA sites DataFrame
        cpa_site = get_cpa_sites(is_polyA, end)

        # Append the read information
        #read_ids.append(read_id)
        reads.append([read_id, umis, cb, chr_name, start, end, strand, cigar, is_polyA, len_pA, int(cpa_site) if cpa_site is not None else None])

    # Create DataFrame from the collected reads
    df = pd.DataFrame(reads, columns=['read_id', 'UMI', 'CB', 'chr', 'start', 'end', 'strand', 'CIGAR', 'is_polyA', 'len_pA', 'cpa_site']) #, index = read_ids)

    # assign reads to read sets
    df = get_readsets(df)

    #df.index = read_ids

    return df

def count_A(sub_sequence):
    """
    Parameters
    ----------
    sub_sequence : string
        a soft clipped region of a read from 5' -> 3' (left to right).
        This is original sequence. (regardless of the read direction)
    
    Returns
    -------        
    number_A : int
        The number of "A"s in the sub_sequence.
    """       
    list_count = [1 if elem == 'A' else 0 for elem in sub_sequence]
    number_A = sum(list_count)
    return number_A


def check_polyA(read, left_end, right_end, rev, percentage_threshold, length_threshold, use_fc):
    """
    Parameters
    ----------
    read : pysam object
        a deduplicated read of interest.

    left_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (-) strand of the genome.
        
        if left_end[0] == 4 -> there is a soft clipped region in the left side of a read.
        left_end[1] -> gives you length of the soft clipped region in the left side of a read.
        
    right_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (+) strand of the genome.   
        
        if right_end[0] == 4 -> there is a soft clipped region in the right side of a read.
        right_end[1] -> gives you length of the soft clipped region in the right side of a read.
        
    percentage_threshold : int
        a percentage threshold for a deduplicated read to be considered as polyA read.
        percentage of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        
    length_threshold : int
        a length threshold for a deduplicated read to be considered as polyA read.
        The number of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        (Note: it does not have to be consecutive number of "A"s)
    
    use_FC : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.        
    
    Returns
    -------        
    True if a read has polyA tail.
    False if a read does not have a polyA tail.
    """          
    full_sequence = read.get_forward_sequence()

    # if read is mapped to negative strand, polyA tail should be on the left end
    if rev == True and left_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] : len(full_sequence)]
            len_pA = left_end[1]
        
        elif use_fc:
            OCS = read.get_tag('XO')
            FCS = read.get_tag('XF')
            difference = OCS - FCS
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] + difference : len(full_sequence)]
            # length of a softclipped region
            len_pA = left_end[1] - difference
            print('potential_polyA: ' + str(potential_polyA))
            print('length: ' + str(len_pA))
            
        # you should use num_A not num_T because you use full_sequence rather than fasta.fetch()
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        
        # for softclipped length of <=6, you want to have 100% "A"        
        if len_pA <= 5:
            percentage_threshold = 100
                
        # decide whether a read is polyA or not
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True, len_pA
        
        else:
            return False, 0
        
    # if right_end[0] == 4, it means soft clipped on the right side of a read.
    # right_end[1] gives how many bases are soft clipped on the right side.
    elif rev == False and right_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] : len(full_sequence)]
            len_pA = right_end[1]
        
        elif use_fc:
            OCS = read.get_tag('XO')
            FCS = read.get_tag('XF')
            difference = FCS - OCS
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] + difference : len(full_sequence)]
            # length of a softclipped region
            len_pA = right_end[1] - difference
            print('potential_polyA: ' + str(potential_polyA))
            print('length: ' + str(len_pA))
            
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        
        # for softclipped length of <=6, you want to have 100% "A"  
        if len_pA <= 5:
            percentage_threshold = 100
        
        # decide whether a read is polyA or not
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True, len_pA
        
        else:
            return False, 0
        
    # a read does not have softclipped region
    # or it has soft clipped region but not at the right direction.
    else:
        return False, 0


def get_cpa_sites(is_polyA, end):    
    """
    Parameters
    ----------
    is_polyA: bool
        Indicates whether the read contains a polyA tail.
    
    start: int
        The start coordinate of the read.
    
    end: int
        The end coordinate of the read.
    
    CIGAR: str
        The CIGAR string representing the alignment.

    df_cpa_sites: pd.DataFrame
        Existing DataFrame to store cpa sites and supporting read counts.

    Returns
    -------            
    df_cpa_sites : pd.DataFrame
        Updated DataFrame with new cpa sites.

    cpa_site : int or None
        the genomic coordinate of the polyA site with the most supporting polyA reads
    """       

    # cpa_site = end of read - softclipped length
    # account for strand direction!!!!!!1
    # check if cpa_site is already in cpa_sites_list
    # if yes, increase supporting_read_count by 1
    # if no, make a new tuple (cpa_site, supporting_read_count) and add it to cpa_sites_list

    if not is_polyA:
        cpa_site = None

    if is_polyA:
        cpa_site = end
        
    return cpa_site


def get_readsets(df):
    """
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing read information including UMI, CB, and is_polyA columns.

    Returns
    -------
    df : pd.DataFrame
        Updated DataFrame with an additional column 'is_polyA_RS' indicating if the read belongs to a polyA read set    
    """

    grouped = df.groupby(['UMI', 'CB'])

    polyA_status = grouped['is_polyA'].any().reset_index()
    polyA_status.rename(columns={'is_polyA': 'is_polyA_RS'}, inplace=True)

    df = df.merge(polyA_status, on=['UMI', 'CB'], how='left')

    return df