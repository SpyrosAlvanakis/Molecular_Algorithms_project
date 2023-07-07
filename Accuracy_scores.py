def check_overlap(range1_start, range1_end, range2_start, range2_end):
    overlap = max(0, min(range1_end, range2_end) - max(range1_start, range2_start) + 1)
    return overlap


def nucleotide_metrics(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true, seq_len):
    true_positives = check_overlap(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true)
    true_negatives = check_overlap(0, start_pos_pred - 1, 0, start_pos_true - 1) + \
                     check_overlap(0, start_pos_pred - 1, end_pos_true + 1, seq_len) + \
                     check_overlap(end_pos_pred + 1, seq_len, 0, start_pos_true - 1) + \
                     check_overlap(end_pos_pred + 1, seq_len, end_pos_true + 1, seq_len)
    false_negatives = check_overlap(0, start_pos_pred - 1, start_pos_true, end_pos_true) + \
                     check_overlap(end_pos_pred + 1, seq_len, start_pos_true, end_pos_true)
    false_positives = check_overlap(start_pos_pred, end_pos_pred, 0, start_pos_true - 1) + \
                      check_overlap(start_pos_pred, end_pos_pred, end_pos_pred + 1, seq_len)

    try:
        performance_coefficient = true_positives / (true_positives + false_positives + false_negatives)
    except ZeroDivisionError:
        performance_coefficient = 0

    try:
        sensitivity = true_positives / (true_positives + false_negatives)
    except ZeroDivisionError:
        sensitivity = 0

    try:
        specificity = true_positives / (true_positives + false_positives)
    except ZeroDivisionError:
        specificity = 0
    return performance_coefficient, sensitivity, specificity


def nucleotide_accuracy(pred_df, mode='nPC', score_group=None):
    if score_group == None:
        pred_df = pred_df.copy()
    else:
        mask = pred_df['Subgroup'].str.contains('_'+str(score_group)+'_')
        pred_df = pred_df[mask].copy()
        pred_df.reset_index(inplace=True, drop=True)

    # filenames = pred_df['Subgroup'].str.split('_').str[0]
    # total_motif_groups_num = len(set(filenames))
    # print(f'Motif groups {total_motif_groups_num}')
    #
    # temp_df_non_zero = pred_df[pred_df[mode] != 0].copy()
    # sequences = pred_df['Subgroup'].str.split('_').str[2]
    # total_sequences_num = len(set(sequences))
    # non_zero_sequence_num = len(set(temp_df_non_zero['Subgroup'].str.split('_').str[2]))
    # print(f'Sequences {total_sequences_num}, {non_zero_sequence_num}')

    total_sites_num = pred_df.shape[0]
    # non_zero_sites_num = temp_df_non_zero.shape[0]
    # print(f'Sites {total_sites_num}')

    if mode in ['nSp', 'nSn', 'nPC']:
        sums = int(pred_df[mode].sum(axis=0))
    else:
        sums = int(pred_df['nPC'].sum(axis=0))

    accuracy = sums / total_sites_num
    return accuracy


def site_metrics():
    pass
