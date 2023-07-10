import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np 
from collections import Counter
from find_file_names import return_the_files
#
#
# def return_the_files(path):
#     txt_files = [file for file in os.listdir(path) if file.endswith(".txt")]
#     only_results = []
#     for file in txt_files:
#         if '_' in file:
#             name_parts = file.split("_")
#             if len(name_parts) > 1:
#                 name = name_parts[0]
#                 only_results.append(name)
#     name_counts = Counter(only_results)
#     unique_names = [name for name, count in name_counts.items() if count > 1]
#     return unique_names


def filter_dataframes(df_list, column_name, value):
    new_dfs = []  # Create an empty list to store the new dataframes
    for df in df_list:
        filtered_df = df[df[column_name] == value]  # Filter the dataframe
        new_dfs.append(filtered_df)  # Append the filtered dataframe to the list
    return new_dfs


def union_of_column_values(df_list, column_name):
    values = set()  # Initialize an empty set to store the unique values
    for df in df_list:
        unique_values = df[column_name].unique()  # Get unique values in column
        values.update(unique_values)  # Add unique values to the set
    return list(values)


def assign_score_group(df, col_name, div_number):
    """This function takes as an input a Dataframe,
    the column of interest, and a number of division.
    It separates the Dataframe in as many parts the number
    of division indicates with respect to the values of the 
    column name and labels them. 
    The group 5 includes the higher values, and the group 1
    the lower values.
    Returns:
        DataFrame: A Dataframe with an extra column that 
        labels the row to a group.
    """
    # Sort dataframe by Score column in descending order
    df = df.sort_values(by=col_name, ascending=False)
    # Reset index after sorting
    df = df.reset_index(drop=True)
    # Use pd.qcut() to create Score_group column
    df['Score_group'] = pd.qcut(df[col_name], q=div_number, labels=range(1, div_number+1))
    return df


def concat_dataframes(df_list):
    result_df = pd.concat(df_list)
    return result_df


def plot_motif_scores(files, id_df_list):
    file_names = [file.split('.')[0] for file in files]
    for df, file_name in zip(id_df_list, file_names):
        plt.figure(figsize=(10, 6))
        for group, group_df in df.groupby("Score_group"):
            plt.hist(group_df["Score"], bins=10, alpha=0.5, label=f"Scoring Group {group}")
        plt.xlabel("Motif Score")
        plt.ylabel("Frequency")
        plt.title(f"Distribution of Motif Scores within Scoring Groups for {file_name}")
        plt.legend()
        plt.show()


def assign_subgroup(df_list, file_name_col, seq_id_col, score_group_col):
    new_df_list = []

    for df in df_list:
        df['Subgroup'] = df[file_name_col].astype(str) + "_" + df[score_group_col].astype(str) + "_" + df[seq_id_col].astype(str)
        new_df_list.append(df)
    return new_df_list


def group_collected_results(results_dir):
    files_raw = [file for file in os.listdir(results_dir) if
                 file.endswith('.csv') and file != 'streme_sites.csv']  # do not include streme results

    dfs_list = []
    for file in files_raw:
        df = pd.read_csv(os.path.join(results_dir, file), index_col=0)
        df = assign_score_group(df, 'Score', 5)
        df['Algorithm'] = file.split('_')[0]  # Add new column with the part before the underscore
        vars()[file.split('.')[0]] = df
        dfs_list.append(vars()[file.split('.')[0]])

    id_df_list = assign_subgroup(dfs_list, 'File_name', 'Sequence_ID', 'Score_group')
    plot_motif_scores(files_raw, id_df_list)

    concatenated_df = pd.concat(id_df_list)
    return concatenated_df