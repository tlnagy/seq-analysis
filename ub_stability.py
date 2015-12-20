
import numpy as np
import scipy.stats
import pandas as pd
import re


def load_stability_data(path='data/ubiquitin_complexes.tsv', mutation_df=None):
    """
    Loads Ubiquitin Rosetta stability calculations as a Pandas dataframe, optionally merging on a per-mutant dataframe.
    :param path: the path to the stability .tsv file
    :param mutation_df: An optional mutation dataframe to merge on; an outer left join will be performed on this dataframe
    :return:
    """
    stability = pd.read_table(path)
    stability = stability[stability['ddg_affinity'].map(lambda x: x.lower().strip() != 'none')]
    stability['ddg_affinity'] = stability['ddg_affinity'].astype(float)
    stability['ddg_stability'] = stability['ddg_stability'].astype(float)
    stability['rank'] = stability['rank'].astype(int)
    pat = re.compile("[0-9]+")
    stability['position'] = stability['mutation'].map(lambda s: int(pat.search(s).group(0)))
    stability = stability.rename(columns = {'position': 'positions', 'mutant_aa': 'amino acids'})
    stability['positions'] = stability['positions'].astype(int)
    stability = stability.set_index(['positions', 'amino acids'])
    if mutation_df is not None:
        stability = stability.combine_first(mutation_df)
    return stability
