import pandas as pd 
import math
import sys


def select_algorithms_to_vote(grouped_df, algs):
    mask_algs = grouped_df['Algorithm'].isin(algs)
    grouped_df = grouped_df[mask_algs].copy()
    grouped_df.reset_index(inplace=True, drop=True)
    return grouped_df


def calculate_voting_results(grouped_df):

    grouped_df['Width'].fillna(grouped_df['Site'].str.len(), inplace=True)
    grouped_df['Position_Range'] = grouped_df.apply(
        lambda row: range(int(row['Starting_position']), int(row['Starting_position'] + row['Width'])), axis=1)
    # Make sure Position_Range contains list of positions, not range objects
    grouped_df['Position_Range'] = grouped_df['Position_Range'].apply(list)

    # Expand each list of positions into new rows, keeping the same index
    positions = grouped_df['Position_Range'].apply(pd.Series).stack()

    # Reset index and rename columns for clarity
    positions = positions.reset_index().rename(columns={'level_0': 'original_index', 0: 'Position'})

    # Add Subgroup information back into new DataFrame
    positions['Subgroup'] = positions['original_index'].apply(lambda idx: grouped_df.loc[idx, 'Subgroup'])

    # Now, for each subgroup, we can count the occurrences of each position
    votes = positions.groupby(['Subgroup', 'Position']).size().reset_index().rename(columns={0: 'Votes'})

    # Pivot to make each Subgroup a separate column
    voting_results = votes.pivot(index='Position', columns='Subgroup', values='Votes').fillna(0)

    voting_results = voting_results.transpose()
    voting_results.columns = voting_results.columns.astype(int)
    
    return voting_results


def apply_smoothing(grouped_df, width_column, voting_df):
    # Calculate rounded width values
    grouped_df['Rounded_Width'] = grouped_df[width_column].apply(lambda x: math.ceil(x / 2))
    
    # Transpose the DataFrame
    transposed_df = voting_df.transpose()
    
    # Apply rolling window operation row-wise
    window_size = int(grouped_df['Rounded_Width'].median())
    smoothed_df = transposed_df.rolling(window=window_size, center=True, min_periods=1).sum()
    
    # Transpose the result back to the original orientation
    smoothed_df = smoothed_df.transpose()
    
    return smoothed_df


def get_max_and_range(smoothed_df):
    # Create a DataFrame with maximum values and their positions (column names)
    max_values_df = smoothed_df.idxmax(axis=1).to_frame('Position')
    max_values_df['Max_Value'] = smoothed_df.max(axis=1)

    # Convert 'Position' column to integer, as it might be in string format
    max_values_df['Position'] = max_values_df['Position'].astype(int)

    # Create a new column with the range
    max_values_df['Range'] = max_values_df['Position'].apply(lambda x: range(x - 7, x + 8))
    
    return max_values_df


def counts_of_common(grouped_df):
    # Create dictionaries to store the counts of sites and starting positions
    site_counts = {}
    starting_position_counts = {}

    # Iterate over each row in the concatenated dataframe
    for index, row in grouped_df.iterrows():
        subgroup = row['Subgroup']
        site = row['Site']
        starting_position = row['Starting_position']
    
        # Update the site count dictionary
        if pd.notnull(site):
            if subgroup not in site_counts:
                site_counts[subgroup] = {}
            if site not in site_counts[subgroup]:
                site_counts[subgroup][site] = 0
            site_counts[subgroup][site] += 1
    
        # Update the starting position count dictionary
        if pd.notnull(starting_position):
            if subgroup not in starting_position_counts:
                starting_position_counts[subgroup] = {}
            if starting_position not in starting_position_counts[subgroup]:
                starting_position_counts[subgroup][starting_position] = 0
            starting_position_counts[subgroup][starting_position] += 1
    return starting_position_counts, site_counts


