# Calculate p-values for the difference between mutant and WT barcodes

import numpy as np
import scipy.stats
import pandas as pd

#################################
# utility functions
#################################


def _select(df, index, value):
    """Convenience method for indexing in multi-index data frames."""
    return df.iloc[df.index.get_level_values(index) == value]


def _weighted_avgs(x):
    return x["slope"] * np.log(x["t0"]) / np.log(x["t0"]).sum()


def _get_barcodes(df, index, mutant_type):
    """Get the barcode slopes corresponding to a mutant index."""
    pos, aa = index
    r = _select(df, 'positions', pos)
    return _select(r, mutant_type, aa)


#################################
# test functions to choose from
#################################


def perm_test(base_slopes, pertubation_slopes, n_perm_iters=5000, stat=np.mean):
    """
    Calculates a p-value using a permutation test for the 2-sided (absolute) difference between the slopes.
    """
    mutant_len = len(pertubation_slopes)
    real_diff = np.abs(stat(pertubation_slopes) - stat(base_slopes))
    shuffled = np.concatenate((pertubation_slopes, base_slopes))
    np.random.shuffle(shuffled)
    b = 0
    for i in range(0, n_perm_iters):
        np.random.shuffle(shuffled)
        if np.abs(stat(shuffled[:mutant_len]) - stat(shuffled[mutant_len:])) > real_diff: b += 1
    return b / n_perm_iters # leave biased; it makes more sense


def welch_test(base_slopes, pertubation_slopes):
    """
    Calculates a p-value using a 2-sided Welch's t-test.
    """
    return scipy.stats.ttest_ind(pertubation_slopes, base_slopes, equal_var=False)[1]


#################################
# mutant vs. WT
#################################

def _test_wt_mutant(df, barcode_df, plain_wt, weighted_wt, mutant_type, weighted=True, method=perm_test):
    """
    Returns a Bonferroni-corrected p-value column.
    """
    group, day, pos, aa = barcode_df.name
    r = _select(df, 'group', group)
    r = _select(r, 'days', day)
    r = _select(r, 'positions', pos)
    r = _select(r, mutant_type, aa)

    if weighted:
        wt_slopes = weighted_wt
        mutant_slopes = _weighted_avgs(r)
    else:
        wt_slopes = plain_wt
        mutant_slopes = r['slope']
    return method(wt_slopes, mutant_slopes) * len(barcode_df)


def calc_wt_mutant_pvals(df, barcode_df, by_codon=False, weighted=True, method=perm_test, pval_col_name='p-value'):
    """
    Calculates a Bonferroni-corrected p-value for the difference between WT and mutant barcode slopes,
    adding a column with these p-values.
    :param df: The input multiindexed AA/codon-level dataframe as defined in seq_analysis.py
    :param barcode_df: The input multiindexed barcode-level dataframe as defined in seq_analysis.py
    :param by_codon: Group by (codon, position) rather than (AA, position)
    :param weighted: Calculate using weighted slopes
    :param method: A function of (WT slopes, mutant slopes) that returns the uncorrected p-value.
    :param pval_col_name: The name of the column to add to the data frame.
    :return: The mutant-level input data frame 'df' with columns pval_col_name and '-log10(' + pval_col_name + ')' added
    """

    mutant_type = "codons" if by_codon else "amino acids"
    wt = _select(df, 'amino acids', 'WT')

    # define outside of lambda so we only do it once
    plain_wt = wt['slope']
    weighted_wt = _weighted_avgs(wt).values

    barcode_df[pval_col_name] = barcode_df.apply(lambda x: _test_wt_mutant(df, x, plain_wt, weighted_wt, mutant_type, weighted=weighted, method=method), axis=1)
    barcode_df['-log10(' + pval_col_name + ')'] = barcode_df[pval_col_name].map(lambda x: -np.log10(x))
    return barcode_df


#################################
# case vs. control
#################################


def calc_case_control_diff(control_df, case_df, by_codons=False, weighted=True, diff_col_name='case-control slope diff', control_suffix=': control', case_suffix=': case'):
    """
    Calculates the difference (case slopes - control slopes).
    :param control_df:
    :param case_df:
    :param by_codon: Group by (codon, position) rather than (AA, position)
    :param weighted: Calculate using weighted slopes
    :param diff_col_name: The name of the column for the difference to add
    :return: A dataframe with columns from control_df and case_df, where each column is given the name original_col_name + case_suffix
    or original_col_name + control_suffix. An additional column with name diff_col_name is added.
    """
    slope_name = 'weighted mean slope' if weighted else 'unweighted mean slope'
    case_df = case_df.copy()
    control_df = control_df.copy()
    case_df = case_df.reset_index('group', drop=True).reset_index('days', drop=True).reset_index(['amino acids', 'positions'])
    control_df = control_df.reset_index('group', drop=True).reset_index('days', drop=True).reset_index(['amino acids', 'positions'])
    merged = pd.merge(control_df, case_df, how='inner', suffixes=[control_suffix, case_suffix], left_on=['amino acids', 'positions'], right_on=['amino acids', 'positions'])
    merged[diff_col_name] = merged[slope_name + case_suffix] - merged[slope_name + control_suffix]
    return merged.set_index(['positions', 'amino acids'])


def _test_case_control(diff, barcode_control, barcode_case, mutant_type, weighted=True, method=perm_test):
    """
    Returns a Bonferroni-corrected p-value column.
    """
    control_r = _get_barcodes(barcode_control, diff.name, mutant_type)
    case_r = _get_barcodes(barcode_case, diff.name, mutant_type)

    if weighted:
        control_slopes = _weighted_avgs(control_r)
        case_slopes = _weighted_avgs(case_r)
    else:
        control_slopes = control_r['slope']
        case_slopes = case_r['slope']
    return method(control_slopes, case_slopes) * len(diff)


def calc_case_control_pvals(case_control_diff_df, barcode_control, barcode_case, by_codon=False, weighted=True, method=perm_test, pval_col_name='d-value'):
    """
    Calculates a Bonferroni-corrected p-value for the difference between case and control barcode slopes.
    For example, call this to calculate p-values for the difference between two doses of a drug.
    :param case_control_diff_df: The data frame output from calc_wt_mutant_pvals (the names of the columns that method added isn't needed)
    :param barcode_control: The input multiindexed barcode-level dataframe for the control as defined in seq_analysis.py
    :param barcode_case: The input multiindexed barcode-level dataframe for the case as defined in seq_analysis.py
    :param by_codon: Group by (codon, position) rather than (AA, position)
    :param weighted: Calculate using weighted slopes
    :param method: A function of (WT slopes, mutant slopes) that returns the uncorrected p-value
    :param pval_col_name: The name of the column to add to the data frame.
    :return: The mutant-level input data frame 'case_control_diff_df' with columns pval_col_name and '-log10(' + pval_col_name + ')' added
    """

    mutant_type = "codons" if by_codon else "amino acids"

    case_control_diff_df[pval_col_name] = case_control_diff_df.apply(lambda diff_row: _test_case_control(diff_row, barcode_control, barcode_case, mutant_type, weighted=weighted, method=method), axis=1)
    case_control_diff_df['-log10(' + pval_col_name + ')'] = case_control_diff_df[pval_col_name].map(lambda x: -np.log10(x))
    return case_control_diff_df

