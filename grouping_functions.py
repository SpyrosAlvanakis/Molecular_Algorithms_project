import pandas as pd
import os
import numpy as np 

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
    the column of out interest and a number of division.
    It separetes the Dataframe in as many parts the number
    of division indicates with respect to the values of the 
    column name and labels them. 
    The group 5 includes the higher values and the group 1
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
