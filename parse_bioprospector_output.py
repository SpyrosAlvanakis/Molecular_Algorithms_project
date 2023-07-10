import re
import os
import pandas as pd
from find_file_names import return_the_files


def get_site_info_BP(site_info, motif_info, file_name):
    data = []
    motif_id = 0  # Initial motif ID
    initial_position = site_info[0][0]  # Initial position
    score_mapping = {motif[0]: motif[3] for motif in motif_info}  # Map motif_id to score

    for site in site_info:
        site_id, site_number, _, starting_point, motif_sequence = site

        if site_number == '1' and site_id == initial_position:
            motif_id += 1

        data.append(
            [site_id, site_number, starting_point, motif_sequence, score_mapping.get(str(motif_id), None), file_name])

    columns = ['Sequence_ID', 'Site_number', 'Starting_Point', 'Motif_Sequence', 'Score', 'File_Name']
    df_info = pd.DataFrame(data, columns=columns)
    return df_info


def process_Bioprospector_output(directory_of_files):
    motif_pattern = r'Motif\s+#(\d+):\s+\((\w+/\w+)\)\n\*+\nWidth \((\d+), \d+\);\s+Gap \[\d+, \d+\];\s+MotifScore (\d+\.\d+);\s+Sites (\d+)'
    site_pattern = r'>(\d+-\d+-\w+)\s+len\s\d+\s+site\s+#(\d+)\s+(\w+)\s+(\d+)\n(\w+)'

    file_names = return_the_files(directory_of_files)
    file_names = sorted(file_names)
    site_dfs = []

    for name in file_names:
        file_paths = [file for file in os.listdir(directory_of_files) if file.startswith(f"{name}_")]

        for file_path in file_paths:
            with open(os.path.join(directory_of_files, file_path), 'r') as file:
                Bioprospector_output = file.read()

            motif_info = re.findall(motif_pattern, Bioprospector_output)
            site_info = re.findall(site_pattern, Bioprospector_output)
            df_site = get_site_info_BP(site_info, motif_info, name)
            site_dfs.append(df_site)

    site_df = pd.concat(site_dfs, ignore_index=True)

    site_df['Site_number'] = pd.to_numeric(site_df['Site_number'])
    site_df['Starting_Point'] = pd.to_numeric(site_df['Starting_Point'])

    # Calculating the width of the motif as the length of the 'Motif_Sequence'
    site_df['Width'] = site_df['Motif_Sequence'].apply(len)

    # Final result - selected columns and rename
    result_df = site_df[['File_Name', 'Sequence_ID', 'Motif_Sequence', 'Score', 'Starting_Point', 'Width']].rename(
        columns={'File_Name': 'File_name', 'Sequence_ID': 'Sequence_ID', 'Motif_Sequence': 'Site', 'Score': 'Score',
                 'Starting_Point': 'Starting_position', 'Width': 'Width'})

    return result_df