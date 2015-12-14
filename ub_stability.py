
import numpy as np
import scipy.stats
import pandas as pd


def load_stability_data(path='data/ubiquitin_complexes.tsv', mutation_df=None, by_codon=False):
    """
    Loads Ubiquitin Rosetta stability calculations as a Pandas dataframe, optionally merging on a per-mutant dataframe.
    :param path:
    :param mutation_df:
    :return:
    """
    mutant_type = "codons" if by_codon else "amino acids"
    stability = pd.read_table(path)
    stability['position'] = stability['mutation'].map(lambda s: s[3:-1])
    stability.set_index(['PDB', 'position', 'wt_aa', 'mutant_aa'], inplace=True)
    if mutation_df is not None:
        stability = mutation_df.reset_index().merge(how='outer', right=stability.reset_index(), on_left=[mutant_type, 'positions'], on_right=['aa', 'position'])
    return stability
