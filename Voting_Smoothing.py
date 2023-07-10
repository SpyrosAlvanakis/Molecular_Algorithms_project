import pandas as pd 
import math
import sys
import os
import matplotlib.pyplot as plt
import numpy as np


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


# def get_max_and_range(smoothed_df):
#     # Create a DataFrame with maximum values and their positions (column names)
#     max_values_df = smoothed_df.idxmax(axis=1).to_frame('Position')
#     max_values_df['Max_Value'] = smoothed_df.max(axis=1)

#     # Convert 'Position' column to integer, as it might be in string format
#     max_values_df['Position'] = max_values_df['Position'].astype(int)

#     # Create a new column with the range
#     max_values_df['Range'] = max_values_df['Position'].apply(lambda x: range(x - 7, x + 8))
    
#     return max_values_df

def get_max_and_range(smoothed_df):
    # Create a DataFrame with maximum values and their positions (column names)
    max_values_df = smoothed_df.idxmax(axis=1).to_frame('Position')
    max_values_df['Max_Value'] = smoothed_df.max(axis=1)

    # Convert 'Position' column to integer, as it might be in string format
    max_values_df['Position'] = max_values_df['Position'].astype(int)

    # Reset index to get regular row indexes
    max_values_df.reset_index(inplace=True)

    # Rename 'index' column to 'Subgroup'
    max_values_df.rename(columns={'index': 'Subgroup'}, inplace=True)

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


def make_motif_group_df(text):
    headers = []
    sequences = []

    # Iterate over the lines
    for line in text:
        line = line.strip()
        if line.startswith('>'):
            # Extract the header (row starting with ">")
            headers.append(line[1:])
        else:
            # Extract the sequence
            sequences.append(line)

    # Create a DataFrame from the lists
    df = pd.DataFrame({'Sequence_ID': headers, 'Sequence': sequences})
    return df


def add_predicted_site_sequence(predictions_df, input_dataset_path):
    input_files = os.listdir(input_dataset_path)

    predictions_df['Predicted_sequence'] = pd.DataFrame(['' for _ in range(predictions_df.shape[0])])

    for motif_group_file in input_files:
        with open(os.path.join(input_dataset_path, motif_group_file), 'r') as file:
            motif_group = file.readlines()

        motif_group_df = make_motif_group_df(motif_group)
        for _, row in motif_group_df.iterrows():
            sequence = str(row['Sequence']).strip()
            seq_id = row['Sequence_ID']
            mask_predicted_seq_id = predictions_df['Subgroup'].str.contains(seq_id)
            seq_index = list(predictions_df[mask_predicted_seq_id].index)
            if len(seq_index) > 0:
                for idx in seq_index:
                    start_idx = predictions_df.iloc[idx]['Position']-7
                    end_idx = predictions_df.iloc[idx]['Position']+7
                    site_sequence = sequence[start_idx:end_idx+1]
                    predictions_df.at[idx, 'Predicted_sequence'] = site_sequence

    return predictions_df

def plot_starting_position_subplots(starting_position_counts, subgroups):
    num_subplots = len(subgroups)
    num_cols = min(num_subplots, 3)
    num_rows = (num_subplots - 1) // num_cols + 1

    if num_subplots == 1:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.bar(*zip(*starting_position_counts[subgroups[0]].items()))
        ax.set_xlabel('Starting Position')
        ax.set_ylabel('Count')
        ax.set_title(f'Starting Position Counts - Subgroup: {subgroups[0]}')
        ax.tick_params(axis='x', rotation=90)
    else:
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 6*num_rows))
        fig.subplots_adjust(hspace=0.4)

        for i, subgroup in enumerate(subgroups):
            counts = starting_position_counts[subgroup]
            positions = list(counts.keys())
            counts = list(counts.values())

            if num_rows > 1:
                ax = axes[i // num_cols, i % num_cols]
            else:
                ax = axes[i % num_cols]

            ax.bar(positions, counts)
            ax.set_xlabel('Starting Position')
            ax.set_ylabel('Count')
            ax.set_title(f'{subgroup}')
            ax.tick_params(axis='x', rotation=90)

    plt.show()
    
def find_similar_motifs(streme_results, ground_truth):
    similar_motifs = {}

    for _, row_gt in ground_truth.iterrows():
        target_motif = row_gt['Predicted_sequence']
        for _, row_sr in streme_results.iterrows():
            site = row_sr['Site']
            differences = sum(s != t for s, t in zip(target_motif, site))
            if differences <= 8:
                similar_motifs[target_motif] = differences
                break

    return similar_motifs

def add_weights_with_similarity(predictions, streme_results, weight_set):
    # Initialize the 'Weights' column to 1
    predictions['Weights'] = 1

    # Define the weights for different levels of similarity
    similarity_weights = weight_set

    # Get the unique file names in streme_results
    unique_file_names = streme_results['File_name'].str.split('_', n=1, expand=True)[0].unique()

    # Iterate over each unique file name
    for file_name in unique_file_names:
        # Filter the streme_results for the current file name
        file_results = streme_results[streme_results['File_name'].str.startswith(file_name)]

        # Filter the ground_truth based on the file name part
        filtered_ground_truth = predictions[predictions['Subgroup'].str.startswith(file_name)]

        # Find similar motifs using the find_similar_motifs function
        similar_motifs = find_similar_motifs(file_results, filtered_ground_truth)

        # Iterate over each row in the filtered ground_truth
        for _, row in filtered_ground_truth.iterrows():
            # Check if the target motif is in the similar motifs
            if row['Predicted_sequence'] in similar_motifs:
                # Get the number of differences
                differences = similar_motifs[row['Predicted_sequence']]

                # Add the corresponding weight based on the number of differences
                weight = similarity_weights.get(differences, 0)
                predictions.loc[_, 'Weights'] += weight

    return predictions
