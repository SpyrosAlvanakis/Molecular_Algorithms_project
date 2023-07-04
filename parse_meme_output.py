import re
import os
import pandas as pd


def parse_meme_motif(motif_info, filename):
    meme_dict = {'File_name': [],
                 'Motif_ID': [],
                 'Width': [],
                 'Score': [],
                 'Sites': [],
                 'con': []
                 }
    for i, match in enumerate(motif_info):
        # print(match)
        meme_dict['File_name'].append(filename)
        meme_dict['Motif_ID'].append(int(match.group('index')))
        meme_dict['Width'].append(int(match.group('width')))
        meme_dict['Score'].append(float(match.group('evalue')))
        meme_dict['Sites'].append(int(match.group('sites')))
        meme_dict['con'].append(match.group('motif'))
    meme_motif_df = pd.DataFrame(meme_dict)
    # meme_motif_df.sort_values(axis=0, by='Score', inplace=True, ascending=False)
    return meme_motif_df


def parse_meme_sites(sites_info, motif_df, filename):
    meme_site_pattern = r'(?P<seq_id>\d+-\d+-(?:forward|reverse))\s*(?P<start_number>\d+)\s+\S+\s+(?P<site>\S+)\s+\S+\s+\S+'
    # meme_site_pattern = r'(?P<site_id>\S+)\s+(?P<start_number>\d+)\s+\S+\s+(?P<motif_sequence>\S+)\s+\S+\s+\S+'
    meme_site_dict = {'Sequence_ID': [],
                      'Site': [],
                      'Starting_Point': [],
                      'Score': [],
                      'Width': [],
                      'File_Name': []
                      }

    for i, site in enumerate(sites_info):
        # print(site[0])
        # print(site[1])
        # sites_per_motif = re.findall(meme_site_pattern, site[1])
        motif_mask = motif_df['Motif_ID'] == int(site[0])
        motif_index = motif_df[motif_mask].index

        sites_per_motif = re.finditer(meme_site_pattern, site[1])
        # print(sites_per_motif)
        for j, match in enumerate(sites_per_motif):
            meme_site_dict['Sequence_ID'].append(match.group('seq_id'))
            meme_site_dict['Score'].append(float(motif_df['Score'].loc[motif_index]))
            meme_site_dict['Starting_Point'].append(int(match.group('start_number')))
            meme_site_dict['Site'].append(match.group('site'))
            meme_site_dict['Width'].append(int(motif_df['Width'].loc[motif_index]))
            meme_site_dict['File_Name'].append(filename)

    meme_site_df = pd.DataFrame(meme_site_dict)
    return meme_site_df


def parse_meme_files():
    meme_dir = os.path.join(os.getcwd(), 'Results/MEME')
    files = os.listdir(meme_dir)
    filtered_files = list(filter(lambda name: name if name.find('10') == -1 else '', files))
    print(len(filtered_files))

    meme_motif_pattern = r'MOTIF\s+(?P<motif>\w+)\s+MEME-(?P<index>\d)\s+width\s+=\s+(?P<width>\d+)\s+sites\s+=\s+(?P<sites>\d+).+E-value\s*=\s*(?P<evalue>\d+(?:\.\d+)?(?:e[+-]?\d+)?)'
    sites_pattern = r'MEME-(\d+) sites sorted by position p-value\n(?:.*\n){3}((?:(?!-+).*\n)*)'

    results_list = []
    for filename in filtered_files:
        with open(os.path.join(meme_dir, filename), 'r') as file:
            meme_output = file.read()

        meme_motif_info = re.finditer(meme_motif_pattern, meme_output)
        meme_sites = re.findall(sites_pattern, meme_output)

        name = filename.split('_')[0]
        meme_motif_df = parse_meme_motif(meme_motif_info, name)
        meme_sites_df = parse_meme_sites(meme_sites, meme_motif_df, name)

        results_list.append(meme_sites_df)

    results_df = pd.concat(results_list, ignore_index=True)
    return results_df


if __name__ == "__main__":
    data = parse_meme_files()
    data.sort_values(by=['File_Name', 'Sequence_ID', 'Score'], inplace=True, ascending=False)
    data.to_csv('meme_sites.csv', encoding='utf-8')
