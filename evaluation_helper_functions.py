import pandas as pd
import os
import re
import seaborn as sns
from matplotlib import pyplot as plt

from Accuracy_scores import nucleotide_metrics, site_accuracy, nucleotide_accuracy, nucleotide_acc_weights
from Voting_Smoothing import make_motif_group_df


def initialize_required_columns(predictions_df):
    predictions_df['Target_Start_pos'] = pd.DataFrame([0 for _ in range(predictions_df.shape[0])])
    predictions_df['Target_Stop_pos'] = pd.DataFrame([0 for _ in range(predictions_df.shape[0])])
    predictions_df['Input_Sequence'] = pd.DataFrame(['' for _ in range(predictions_df.shape[0])])
    return predictions_df


def add_ground_truth_columns(predictions_df, ground_truth_dir, input_datasets_dir):
    input_files = os.listdir(input_datasets_dir)
    # ground_truth_motif_files = os.listdir(ground_truth_dir)

    predictions_df = initialize_required_columns(predictions_df.copy())
    for motif_group_file in input_files:
        filename = motif_group_file.split('.txt')[0]

        with open(os.path.join(input_datasets_dir, motif_group_file), 'r') as file:
            motif_group = file.readlines()

        motif_group_df = make_motif_group_df(motif_group)

        # read ground truth motif files as csv with custom named columns
        col_names = ['ECK1', 'File_name', 'ECK2', 'Starting_position', 'Ending_position', 'Orientation',
                     'Association', 'Sequence', 'Additional', 'Label']
        targets_df = pd.read_csv(os.path.join(ground_truth_dir, motif_group_file), sep='\t', names=col_names)
        targets_df['Sequence_ID'] = targets_df['Starting_position'].astype(str) + '-' + targets_df[
            'Ending_position'].astype(str) + '-' + targets_df['Orientation']

        for idx, target_motif in targets_df.iterrows():
            # Match ground truth sequence id with the sequence id in the motif group input file to get the original
            # sequence
            seq_id = str(target_motif['Sequence_ID']).strip()
            seq_id_mask = motif_group_df['Sequence_ID'] == seq_id
            seq_index = list(motif_group_df[seq_id_mask].index)
            if len(seq_index) > 0:
                index = seq_index[0]
                input_sequence = str(motif_group_df.iloc[index]['Sequence']).strip()

                target_motif_sequence = str(target_motif['Sequence']).strip()
                uppercase_part = re.findall('[A-Z]+', target_motif_sequence)

                try:
                    target_site = uppercase_part[0]
                    # find the position of the ground truth motif inside the input sequence
                    target_site_start_pos = input_sequence.index(target_site)
                    target_site_end_pos = target_site_start_pos + len(target_site)

                    # Use a mask to match the current sequence id with the sequence ids in the predictions dataframe
                    mask_predicted_seq_id = predictions_df['Subgroup'].str.contains(seq_id)
                    predicted_sites = predictions_df[mask_predicted_seq_id]
                    for pred_idx, row in predicted_sites.iterrows():
                        # Store in the predictions matched, the target positions as found above and the original input
                        # sequence
                        predictions_df.at[pred_idx, 'Target_Start_pos'] = target_site_start_pos
                        predictions_df.at[pred_idx, 'Target_Stop_pos'] = target_site_end_pos
                        predictions_df['Input_Sequence'] = input_sequence
                except IndexError:
                    print(filename, seq_id)
            else:
                pass

    predictions_df['Predicted_Start_pos'] = predictions_df['Position'] - 7
    predictions_df['Predicted_Stop_pos'] = predictions_df['Position'] + 7
    return predictions_df


def add_nucleotide_metrics_columns(predictions_df):
    predictions_df['nPC'] = pd.DataFrame([0 for _ in range(predictions_df.shape[0])])
    predictions_df['nSn'] = pd.DataFrame([0 for _ in range(predictions_df.shape[0])])
    predictions_df['nSp'] = pd.DataFrame([0 for _ in range(predictions_df.shape[0])])

    for i, row in predictions_df.iterrows():
        start_pred = row['Predicted_Start_pos']
        stop_pred = row['Predicted_Stop_pos']
        start_true = row['Target_Start_pos']
        stop_true = row['Target_Stop_pos']
        seq_len = len(row['Input_Sequence']) - 1
        npc, nsn, nsp = nucleotide_metrics(start_pred, stop_pred, start_true, stop_true, seq_len)
        predictions_df.at[i, 'nPC'] = npc
        predictions_df.at[i, 'nSn'] = nsn
        predictions_df.at[i, 'nSp'] = nsp

    return predictions_df


def fix_alg_names(algs):
    alg_names = {'bioprospector': 'BP',
                 'meme': 'ME',
                 'motifSampler': 'MS',
                 'mdscan': 'MD'
                 }

    alg_abbr = [alg_names[name] for name in algs]
    names_string = '-'.join(alg_abbr)
    return names_string


def plot_nucleotide_accuracies(predictions, mode='nPC', dataset='Type B', margin='200', algorithms=None, runs='10'):
    if algorithms is None:
        algorithms = ['bioprospector', 'meme', 'motifSampler', 'mdscan']
    algorithms = fix_alg_names(algorithms)
    total = round(nucleotide_accuracy(predictions, mode=mode), 3)
    group1 = round(nucleotide_accuracy(predictions, mode=mode, score_group=1), 3)
    group2 = round(nucleotide_accuracy(predictions, mode=mode, score_group=2), 3)
    group3 = round(nucleotide_accuracy(predictions, mode=mode, score_group=3), 3)
    group4 = round(nucleotide_accuracy(predictions, mode=mode, score_group=4), 3)
    group5 = round(nucleotide_accuracy(predictions, mode=mode, score_group=5), 3)

    results = {'Total': total,
               'Score Group 1': group1,
               'Score Group 2': group2,
               'Score Group 3': group3,
               'Score Group 4': group4,
               'Score Group 5': group5
               }

    plt.rcParams["figure.autolayout"] = True
    plt.figure(figsize=(10,8))
    sns.set_style('whitegrid')
    ax = sns.barplot(x=list(results.keys()), y=list(results.values()), color='blue', palette='hls')
    ax.bar_label(ax.containers[0], fmt='%.3f')

    plt.xlabel('Score Groups', fontsize=15)
    plt.ylabel(f'Accuracy {mode}', fontsize=15)


    plt.suptitle(f'EMD {algorithms} Accuracy {mode} on Dataset {dataset} with Margin {margin} and {runs} Runs', fontsize=15)
    plt.savefig(f'accuracy_{mode}_{algorithms}_{dataset}{margin}_{runs}.png')
    plt.show()


def plot_site_accuracies(predictions, mode='sPC', dataset='Type B', margin='200', algorithms=None, runs='10'):
    if algorithms is None:
        algorithms = ['bioprospector', 'meme', 'motifSampler', 'mdscan']
    algorithms = fix_alg_names(algorithms)
    total = round(site_accuracy(predictions, mode=mode), 3)
    group1 = round(site_accuracy(predictions, mode=mode, score_group=1), 3)
    group2 = round(site_accuracy(predictions, mode=mode, score_group=2), 3)
    group3 = round(site_accuracy(predictions, mode=mode, score_group=3), 3)
    group4 = round(site_accuracy(predictions, mode=mode, score_group=4), 3)
    group5 = round(site_accuracy(predictions, mode=mode, score_group=5), 3)

    results = {'Total': total,
               'Score Group 1': group1,
               'Score Group 2': group2,
               'Score Group 3': group3,
               'Score Group 4': group4,
               'Score Group 5': group5
               }

    plt.rcParams["figure.autolayout"] = True
    plt.figure(figsize=(10, 8))
    sns.set_style('whitegrid')
    ax = sns.barplot(x=list(results.keys()), y=list(results.values()), color='green', palette='hls')
    ax.bar_label(ax.containers[0], fmt='%.3f')

    plt.xlabel('Score Groups', fontsize=15)
    plt.ylabel(f'Accuracy {mode}', fontsize=15)


    plt.suptitle(f'EMD {algorithms} Accuracy {mode} on Dataset {dataset} with Margin {margin} and {runs} Runs', fontsize=15)
    plt.savefig(f'accuracy_{mode}_{algorithms}_{dataset}{margin}_{runs}.png')
    plt.show()


def plot_nucleotide_weighted_accuracies(predictions, wset=1, mode='nPC', dataset='Type B', margin='200', algorithms='BP-MD-ME-MS', runs='10'):
    total = round(nucleotide_acc_weights(predictions, mode=mode), 3)
    group1 = round(nucleotide_acc_weights(predictions, mode=mode, score_group=1), 3)
    group2 = round(nucleotide_acc_weights(predictions, mode=mode, score_group=2), 3)
    group3 = round(nucleotide_acc_weights(predictions, mode=mode, score_group=3), 3)
    group4 = round(nucleotide_acc_weights(predictions, mode=mode, score_group=4), 3)
    group5 = round(nucleotide_acc_weights(predictions, mode=mode, score_group=5), 3)

    results = {'Total': total}

    if group1 != 0:
        results['Score Group 1'] = group1
    if group2 != 0:
        results['Score Group 2'] = group2
    if group3 != 0:
        results['Score Group 3'] = group3
    if group4 != 0:
        results['Score Group 4'] = group4
    if group5 != 0:
        results['Score Group 5'] = group5

    plt.rcParams["figure.autolayout"] = True
    plt.figure(figsize=(10,8))
    sns.set_style('whitegrid')
    ax = sns.barplot(x=list(results.keys()), y=list(results.values()), color='blue', palette='hls')
    ax.bar_label(ax.containers[0], fmt='%.3f')

    plt.xlabel('Score Groups', fontsize=15)
    plt.ylabel(f'Accuracy {mode}', fontsize=15)

    plt.suptitle(f'EMD {algorithms} Accuracy {mode} on Weighted Dataset {dataset} (weight set {wset}) with Margin {margin} and {runs} Runs', fontsize=15)
    plt.savefig(f'accuracy{mode}{algorithms}{dataset}{margin}_{runs}_wset{wset}.png')
    plt.show()