from __future__ import print_function, division
import warnings
import sys, os
try:
    import platform
    if platform.python_version().startswith('2'):
        warnings.warn("Python 2 is old and yucky, please use 3", RuntimeWarning)
except ImportError:
    sys.exit(2)
import pandas as pd
import numpy as np


def calc_pvals(cc, k = 800):
    cc['p-value'] = np.zeros(len(cc))
    control = cc['slope_control'].values
    case = cc['slope_case'].values
    for i in range(0, k):
        np.random.shuffle(control)
        np.random.shuffle(case)
        rev = 1 if np.random.rand() < 0.5 else -1
        cc['p-value'] = cc['p-value'] + (rev * (case - control) > cc['slope_diff']).astype(int)
    cc['p-value'] = cc['p-value'] / k
    return cc


def calc_diff(fitnesses, group, control_day = 'd1', case_day='d2', n_sims = 800):
    '''
    Calculates the difference between case and control fitness values with p-values calculated by permuting the data across both.
    :param fitnesses: Pandas dataframe from seq_analysis of fitness values at either the codon or amino acid level
    :param group: The name of the group (e.g. PyND); set to None to calculate for all groups
    :param control_day: The name of the day corresponding to the control
    :param case_day: The name of the day corresponding to the perturbation
    :n_simulations: The number of iterations to use in the permutation test; set to 0 to disable
    :return: A Pandas dataframe with columns [slope_control, slope_case, slope_diff, std_control, std_case, p_val]
    '''
    if group is not None:
        days = fitnesses[fitnesses.index.get_level_values(level = 'group') == group]
    else:
        days = fitnesses
    days = days.groupby(level = 'days')
    control = days.get_group(control_day)
    case = days.get_group(case_day)
    case.columns = case.columns.get_level_values(1)
    control.columns = control.columns.get_level_values(1)
    control = control.reset_index()
    case = case.reset_index()
    wanted_cols = [c for c in ['codons', 'positions', 'amino acids', '# unique barcodes', 'weighted mean slope', 'stddev of slope', 'sum of t0 reads'] if c in case.columns]
    wanted_indices = [c for c in ['codons', 'positions', 'amino acids'] if c in case.columns]
    case = case[wanted_cols]
    control = control[wanted_cols]
    merged = pd.merge(control, case, how='inner', on=wanted_indices, suffixes=['_control', '_case'])
    merged.set_index(['codons', 'amino acids', 'positions'], inplace=True)
    merged['slope_diff'] = merged['slope_case'] - merged['slope_control']
    if n_sims > 0:
        calc_pvals(merged, n_sims)
    return merged
