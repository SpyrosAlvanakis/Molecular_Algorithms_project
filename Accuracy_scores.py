def check_overlap(range1_start, range1_end, range2_start, range2_end):
    overlap = max(0, min(range1_end, range2_end) - max(range1_start, range2_start) + 1)
    return overlap


def nucleotide_metrics(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true, seq_len):
    true_positives = check_overlap(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true)
    # true_negatives = check_overlap(0, start_pos_pred - 1, 0, start_pos_true - 1) + \
    #                  check_overlap(0, start_pos_pred - 1, end_pos_true + 1, seq_len) + \
    #                  check_overlap(end_pos_pred + 1, seq_len, 0, start_pos_true - 1) + \
    #                  check_overlap(end_pos_pred + 1, seq_len, end_pos_true + 1, seq_len)
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
        mask = pred_df['Subgroup'].str.contains('_' + str(score_group) + '_')
        pred_df = pred_df[mask].copy()
        pred_df.reset_index(inplace=True, drop=True)

    total_sites_num = pred_df.shape[0]

    if mode in ['nSp', 'nSn', 'nPC']:
        sums = int(pred_df[mode].sum(axis=0))
    else:
        sums = int(pred_df['nPC'].sum(axis=0))

    accuracy = sums / total_sites_num
    return accuracy


def site_accuracy(pred_df, threshold=1, mode='sPC', score_group=None):
    if score_group == None:
        pred_df = pred_df.copy()
    else:
        mask = pred_df['Subgroup'].str.contains('_' + str(score_group) + '_')
        pred_df = pred_df[mask].copy()
        pred_df.reset_index(inplace=True, drop=True)

    true_positives = 0
    false_positives = 0

    target_sequence_overlap = {}
    for idx, row in pred_df.iterrows():
        start_pred = row['Predicted_Start_pos']
        stop_pred = row['Predicted_Stop_pos']
        start_true = row['Target_Start_pos']
        stop_true = row['Target_Stop_pos']

        pred_width = stop_pred-start_pred+1
        seq_id = row['Subgroup'].split('_')[2]

        overlap = check_overlap(start_pred, stop_pred, start_true, stop_true)
        if pred_width - overlap <= threshold:
            true_positives += 1
            if seq_id in list(target_sequence_overlap.keys()):
                target_sequence_overlap[seq_id] += 1
            else:
                target_sequence_overlap[seq_id] = 1
        else:
            false_positives += 1
            if seq_id not in list(target_sequence_overlap.keys()):
                target_sequence_overlap[seq_id] = 0

    false_negatives = list(target_sequence_overlap.values()).count(0)

    try:
        sPC = true_positives / (true_positives + false_positives + false_negatives)
    except ZeroDivisionError:
        sPC = 0

    try:
        sSn = true_positives / (true_positives + false_negatives)
    except ZeroDivisionError:
        sSn = 0

    try:
        sSp = true_positives / (true_positives + false_positives)
    except ZeroDivisionError:
        sSp = 0

    total_motif_groups_num = len(set(pred_df['Subgroup'].str.split('_')[0]))
    if mode == 'sSp':
        accuracy = sSp / total_motif_groups_num
    elif mode == 'sSn':
        accuracy = sSn / total_motif_groups_num
    else:
        accuracy = sPC / total_motif_groups_num

    return accuracy

