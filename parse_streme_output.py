import os
import re
import pandas as pd
import itertools
from find_file_names import return_the_files


def get_best_score_streme(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('letter-probability matrix'):
                parts = line.split()
                for i, part in enumerate(parts):
                    if 'S=' in part:
                        try:
                            e_value = float(parts[i+1])
                            return e_value
                        except ValueError:
                            print(f"Cannot convert '{parts[i+1]}' to float")
                            return None
    return None

# Dictionary to map IUPAC nucleotide codes to the possible bases they represent
IUPAC_dict = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


def generate_sequences_streme(motif):
    # Generate a list of lists where each inner list contains the possible bases for each position in the motif
    bases = [IUPAC_dict[char] for char in motif]
    # Generate all combinations of the possible bases
    combinations = list(itertools.product(*bases))
    # Join each combination into a string and return a list of all possible sequences
    sequences = [''.join(combination) for combination in combinations]
    return sequences


def generate_motif_sequences(file_dir):
    # Get the list of file names from the given directory
    file_names = return_the_files(file_dir)
    # Sort the file names for consistency
    file_names = sorted(file_names)
    # Initialize the list that will hold all results
    all_results = []

    # Loop over each file name
    for file_name in file_names:
        # Loop over each file in the directory
        for file in os.listdir(file_dir):
            # If the file starts with the current file name
            if file.startswith(file_name):
                # Create the full file path
                file_path = os.path.join(file_dir, file)
                # Get the best score from the file
                best_score = get_best_score_streme(file_path)
                # Open the file
                with open(file_path, 'r') as open_file:
                    # Read the contents of the file
                    output_streme = open_file.read()
                # Define the regex pattern to extract motifs
                motif_pattern = r"MOTIF (\d+)-([A-Z]+) STREME-\d+"
                # Find all matches in the file contents
                motifs = re.findall(motif_pattern, output_streme)
                # Create a DataFrame from the motifs
                temp_motifs = pd.DataFrame(motifs, columns=['Motif_ID', 'Con'])
                # Convert the 'Motif_ID' column to integer
                temp_motifs['Motif_ID'] = temp_motifs['Motif_ID'].astype(int)
                # Extract the first part of the file name (before '_') and assign it to the 'File_Name' column
                temp_motifs['File_name'] = file.split('_')[0]
                # Assign the best score to the 'Best_Score' column
                temp_motifs['Score'] = best_score

                # Loop over each row in the DataFrame
                for _, row in temp_motifs.iterrows():
                    # Generate all possible sequences from the 'Con' column
                    sequences = generate_sequences_streme(row['Con'])
                    # For each sequence, append a new row to the results
                    for sequence in sequences:
                        all_results.append({'File_name': file.split('_')[0], 'Site': sequence, 'Score': best_score,
                                            'Width': len(sequence)})

    streme_motifs = pd.DataFrame(all_results)
    return streme_motifs


