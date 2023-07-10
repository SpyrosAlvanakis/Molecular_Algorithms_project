import os
from collections import Counter


def return_the_files(path):
    txt_files = [file for file in os.listdir(path) if file.endswith(".txt")]
    only_results = []
    for file in txt_files:
        if '_' in file:
            name_parts = file.split("_")
            if len(name_parts) > 1:
                name = name_parts[0]
                only_results.append(name)
    name_counts = Counter(only_results)
    unique_names = [name for name, count in name_counts.items() if count > 1] 
    return unique_names