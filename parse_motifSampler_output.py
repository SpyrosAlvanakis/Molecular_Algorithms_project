import re
import os
import pandas as pd


def parse_motif(motif_info, filename):
    motif_dict = {'motif_id': [],
                  'score': [],
                  'filename': []
                  }

    for match in motif_info:
        id = match.group('motif_id')
        score = float(match.group('score'))

        motif_dict['motif_id'].append(id)
        motif_dict['score'].append(score)
        motif_dict['filename'].append(filename)

    motif_df = pd.DataFrame(motif_dict)
    return motif_df


def keep_top5_motifs(df):
    df = df.sort_values(by='score', ascending=False)
    df.reset_index(inplace=True, drop=True)
    new_df = df.loc[:4]
    # print(new_df.shape)
    return new_df


def parse_sites(sites_info, motif_df, filename):
    ms_site_dict = {'Sequence_ID': [],
                    'Site': [],
                    'Starting_position': [],
                    'Score': [],
                    'Width': [],
                    'File_name': []
                    }

    valid_motif_ids = motif_df['motif_id'].tolist()
    for match_sites in sites_info:
        seq_id = match_sites.group('seq_id')
        start_pos = int(match_sites.group('start_pos'))
        motif_id = match_sites.group('motif_id')
        site = match_sites.group('site')
        width = len(site)

        if motif_id in valid_motif_ids:
            motif_mask = motif_df['motif_id'] == motif_id
            motif_index = motif_df[motif_mask].index
            score = float(motif_df['score'].loc[motif_index])

            ms_site_dict['Sequence_ID'].append(seq_id)
            ms_site_dict['Site'].append(site)
            ms_site_dict['Starting_position'].append(start_pos)
            ms_site_dict['Score'].append(score)
            ms_site_dict['Width'].append(width)
            ms_site_dict['File_name'].append(filename)

    sites_df = pd.DataFrame(ms_site_dict)
    return sites_df


def parse_motifSampler_files(ms_dir):
    files = os.listdir(ms_dir)
    filtered_files = list(filter(lambda name: name if name.find('.txt') > -1 else '', files))
    # print(len(filtered_files))

    pattern_id_cs = r'#id: (?P<motif_id>\S+).*cs: (?P<score>\S+)'
    pattern_sites = r'(?P<seq_id>\d+-\d+-(?:forward|reverse)).*misc_feature\s*(?P<start_pos>\d+).*id "(?P<motif_id>\S+)"; site "(?P<site>\S+)";'

    results_list = []
    for filename in filtered_files:
        with open(os.path.join(ms_dir, filename), 'r') as file:
            ms_output = file.read()

        motif_info = re.finditer(pattern_id_cs, ms_output)
        sites_info = re.finditer(pattern_sites, ms_output)

        name = filename.split('_')[0]
        motif_df = parse_motif(motif_info, name)
        motif_df = keep_top5_motifs(motif_df.copy())

        sites_df = parse_sites(sites_info, motif_df, name)
        results_list.append(sites_df)

    results_df = pd.concat(results_list, ignore_index=True)
    # print(results_df.head())
    return results_df


if __name__ == "__main__":
    ms_dir = os.path.join(os.getcwd(), 'Results/MotifSampler')
    data = parse_motifSampler_files(ms_dir)

    data.sort_values(by=['File_name', 'Sequence_ID', 'Score'], inplace=True, ascending=False)
    data.reset_index(drop=True, inplace=True)
    data.to_csv('motifSampler_sites.csv', encoding='utf-8')