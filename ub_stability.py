
import numpy as np
import scipy.stats
import pandas as pd


def load_stability_data(path='data/ubiquitin_complexes.tsv', mutation_df=None):
    """
    Loads Ubiquitin Rosetta stability calculations as a Pandas dataframe, optionally merging on a per-mutant dataframe.
    :param path:
    :param mutation_df:
    :return:
    """
    stability = pd.read_table(path)
    stability['position'] = stability['mutation'].map(lambda s: s[3:-1])
    stability.set_index(['PDB', 'position', 'wt_aa', 'mutant_aa'], inplace=True)
    if mutation_df is not None:
        stability = mutation_df.reset_index().merge(how='right', right=stability.reset_index(), left_on=["amino acids", 'positions'], right_on=['mutant_aa', 'position'])
        stability = stability.drop('amino acids', axis=1).drop('positions', axis=1)
    return stability
