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


def calc_diff(fitnesses, case_group, control_group="Control", control_day='d1', case_day='d1', n_sims=800, by_codons=False):
    '''
    Calculates the difference between case and control fitness values with p-values calculated by permuting the data across both.
    :param fitnesses: Pandas dataframe from seq_analysis of fitness values at either the codon or amino acid level
    :param group: The name of the group (e.g. PyND); set to None to calculate for all groups
    :param control_day: The name of the day corresponding to the control
    :param case_day: The name of the day corresponding to the perturbation
    :n_simulations: The number of iterations to use in the permutation test; set to 0 to disable
    :return: A Pandas dataframe with columns [slope_control, slope_case, slope_diff, std_control, std_case, p_val]
    '''
    control = fitnesses
    case = fitnesses
    if control_group is not None:
        control = fitnesses[fitnesses.index.get_level_values(level = 'group') == control_group]
    if case_group is not None:
        case = fitnesses[fitnesses.index.get_level_values(level = 'group') == case_group]
    if control_day is not None:
        control = control[fitnesses.index.get_level_values(level = 'days') == control_day]
    if case_day is not None:
        case = control[fitnesses.index.get_level_values(level = 'days') == case_day]
    case.columns = case.columns.get_level_values(1) # flatten
    control.columns = control.columns.get_level_values(1) # flatten
    control = control.reset_index()
    case = case.reset_index()
    wanted_cols = [c for c in ['group', 'days', 'positions', '# unique barcodes', 'weighted mean slope', 'stddev of slope', 'weighted stddev', 'sum of t0 reads', 'codons' if by_codons else 'amino acids'] if c in case.columns]
    wanted_indices = [c for c in ['group', 'days', 'positions', 'codons' if by_codons else 'amino acids'] if c in case.columns]
    case = case[wanted_cols]
    control = control[wanted_cols]
    merged = pd.merge(control, case, how='inner', on=wanted_indices, suffixes=['_control', '_case'])
    merged.set_index(wanted_indices, inplace=True)
    merged['slope_diff'] = merged['slope_case'] - merged['slope_control']
    if n_sims > 0:
        calc_pvals(merged, n_sims)
    return merged
