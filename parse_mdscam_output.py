import re
import os
import pandas as pd
from find_file_names import return_the_files


def get_site_info_MD(site_info, motif_info, file_name):
    data = []
    motif_id = 0  # Initial motif ID
    initial_position = site_info[0][0]  # Initial position
    score_mapping = {motif[0]: motif[2] for motif in motif_info}  # Map motif_id to score

    for site in site_info:
        site_id, site_number, _, starting_point, motif_sequence = site
        width = len(motif_sequence)

        if site_number == '1' and site_id == initial_position:
            motif_id += 1

        data.append([site_id, site_number, starting_point, motif_sequence, score_mapping.get(str(motif_id), None), file_name, width])

    columns = ['Sequence_ID', 'Site_number', 'Starting_Point', 'Motif_Sequence', 'Score', 'File_Name', 'Width']
    df_info = pd.DataFrame(data, columns=columns)
    return df_info


def process_MDScan_output(directory_of_files):
    """
    This function read the txt files found in the given directory and parses them assuming they follow the output format
    of MDscan.

    :param directory_of_files: A string that should be the directory containing the output results of MDScan
    :return:
    """
    # Extract the motifs using regex
    motif_pattern = r"Motif\s+(\d+):\s+Wid\s+(\d+);\s+Score\s+([\d.]+);\s+Sites\s+(\d+);\s+Con\s+([ACGT]+);\s+RCon\s+([ACGT]+)"
    # motif_info = re.findall(motif_pattern, mdscan_output)

    # Extract the site information using regex
    site_pattern = r">(\d+-\d+-(?:forward|reverse))\s+Len\s+\d+\s+Site\s+#(\d+)\s+([fr])\s+(\d+)\n([ACGT]+)"

    file_names = return_the_files(directory_of_files)
    file_names = sorted(file_names)
    site_dfs = []

    for name in file_names:
        file_paths = [file for file in os.listdir(directory_of_files) if file.startswith(f"{name}_")]

        for file_path in file_paths:
            with open(os.path.join(directory_of_files, file_path), 'r') as file:
                mdscan_output = file.read()

            motif_info = re.findall(motif_pattern, mdscan_output)
            site_info = re.findall(site_pattern, mdscan_output)
            df_site = get_site_info_MD(site_info, motif_info, name)
            site_dfs.append(df_site)

    site_df = pd.concat(site_dfs, ignore_index=True)

    site_df['Site_number'] = pd.to_numeric(site_df['Site_number'])
    site_df['Starting_Point'] = pd.to_numeric(site_df['Starting_Point'])

    # Final result - selected columns and rename
    result_df = site_df[['File_Name', 'Sequence_ID', 'Motif_Sequence', 'Score', 'Starting_Point', 'Width']].rename(
        columns={'File_Name': 'File_name', 'Sequence_ID': 'Sequence_ID', 'Motif_Sequence': 'Site', 'Score': 'Score',
                 'Starting_Point': 'Starting_position', 'Width': 'Width'})

    return result_df



