def check_overlap(range1_start, range1_end, range2_start, range2_end):
    overlap = max(0, min(range1_end, range2_end) - max(range1_start, range2_start) + 1)
    return overlap


def nucleotide_metrics(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true, seq_len):
    true_positives = check_overlap(start_pos_pred, end_pos_pred, start_pos_true, end_pos_true)
    true_negatives = check_overlap(0, start_pos_pred-1, 0, start_pos_true-1) + \
                     check_overlap(0, start_pos_pred-1, end_pos_true+1, seq_len) + \
                     check_overlap(end_pos_pred+1, seq_len, 0, start_pos_true - 1) + \
                     check_overlap(end_pos_pred+1, seq_len, end_pos_true+1, seq_len)
    false_negative = check_overlap(0, start_pos_pred-1, start_pos_true, end_pos_true) + \
                     check_overlap(end_pos_pred+1, seq_len, start_pos_true, end_pos_true)
    false_positives = check_overlap(start_pos_pred, end_pos_pred, 0, start_pos_true-1) + \
                      check_overlap(start_pos_pred, end_pos_pred, end_pos_pred+1, seq_len)

    try:
        performance_coefficient = true_positives/(true_positives + false_positives + false_negative)
    except ZeroDivisionError:
        performance_coefficient = 0

    try:
        sensitivity = true_positives/(true_positives + false_negative)
    except ZeroDivisionError:
        sensitivity = 0

    try:
        specificity = true_positives/(true_positives + false_positives)
    except ZeroDivisionError:
        specificity = 0
    return performance_coefficient, sensitivity, specificity


def nucleotide_accuracy():
    pass


def site_metrics():
    pass


