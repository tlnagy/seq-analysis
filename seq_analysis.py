from __future__ import print_function, division
import warnings
import sys, os
try:
    import platform
    if platform.python_version().startswith('2'):
        warnings.warn("Python 2 is old and yucky, please use 3", RuntimeWarning)
except ImportError:
    sys.exit(2)
import difflib
from Bio.Seq import Seq
import pandas as pd
import pickle
import numpy as np
import scipy.stats
import utils


def map_dataset(hdf5_datastorepath, group_name, allele_pkl_path=None):
    """
    Loads raw demultiplexed read h5 file and maps it onto the barcode-mutant
    pairing. Does not process the data further.
    """
    print("Loading dataset...", end="", flush=True)
    raw_barcode_data = pd.read_hdf(hdf5_datastorepath, key="grouped_data")
    fuzzy_matching = {group.lower(): group for group in raw_barcode_data.index.levels[0]}
    closest_matches = difflib.get_close_matches(group_name.lower(), fuzzy_matching.keys())
    try:
        matched_group = fuzzy_matching[closest_matches[0]]
    except IndexError:
        print("Use one of the following: {}".format(list(fuzzy_matching.values())))
        return None, None
    idx = pd.IndexSlice
    raw_barcode_data = raw_barcode_data.loc[idx[matched_group], :].reset_index()
    
    num_unique_reads = len(raw_barcode_data)
    unique_read_count_total = raw_barcode_data["counts"].sum()

    if allele_pkl_path is None:
        allele_pkl_path = os.path.join(*[os.path.split(hdf5_datastorepath)[0], "allele_dic_with_WT.pkl"])
    with open(allele_pkl_path, "rb") as f:
        barcode_mutant_map = pd.DataFrame(pickle.load(f)).T.reset_index()
    print("Done.\n\nMapping...", flush=True)

    barcode_mutant_map.columns = ["barcodes", "positions", "codons"]
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].astype(np.str)

    barcode_mutant_map["WT"] = barcode_mutant_map["codons"] == "WT"
    # add dummy value for WT barcodes
    barcode_mutant_map["amino acids"] = "@"
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).transcribe()))
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "amino acids"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).translate()))

    mapped_barcode_data = raw_barcode_data.merge(barcode_mutant_map, on="barcodes", how="inner")
    mapped_barcode_data["positions"] = mapped_barcode_data["positions"].astype(np.int)
    mapped_barcode_data.index.name = "idxs"

    print("Percent of unique barcodes mapped: {:.2%}".format(len(mapped_barcode_data)/num_unique_reads))
    print("Percent of total reads mapped: {:.2%}\n".format(mapped_barcode_data["counts"].sum()/unique_read_count_total))

    print("group   % of total barcodes found")
    for name, group in mapped_barcode_data.groupby(["days", "timepoints"]):
            print("{}    {:.2%}".format(name, len(group)/len(barcode_mutant_map)))

    return raw_barcode_data, mapped_barcode_data


def process_data(mapped_barcode_data):
    print("\nProcessing data...", flush=True)
    # transform data
    processed_barcodes = mapped_barcode_data.pivot_table(index=["days", "timepoints", "barcodes", "codons","amino acids", "positions"],
                                                         values=["counts", "rel_freq", "rel_wt"])

    processed_barcodes = processed_barcodes.unstack("days").unstack("timepoints")

    idx = pd.IndexSlice
    # Throw out barcodes on a per day basis if they have zero counts either
    before = len(processed_barcodes)

    # consider days independently
    processed_barcodes = processed_barcodes.stack("days")
    processed_barcodes = processed_barcodes[pd.notnull(processed_barcodes.loc[:, idx["counts"]]).sum(axis=1) == len(processed_barcodes.columns)]
    processed_barcodes = processed_barcodes.stack("timepoints").unstack("days").unstack("timepoints")
    print("\nDiscarding {} barcodes on days where they had a count of 0".format(before-len(processed_barcodes)))

    before = len(processed_barcodes)
    # Throw out values that have a count 1 in any experiment on a per day basis
    processed_barcodes = processed_barcodes.stack("days")
    processed_barcodes = processed_barcodes[processed_barcodes[processed_barcodes == 1].count(axis=1) == 0]
    processed_barcodes = processed_barcodes.stack("timepoints").unstack("days").unstack("timepoints")
    print("Discarding {} barcodes on days that had a count of 1 at any timepoint".format(before-len(processed_barcodes)))

    sums = processed_barcodes["counts"].sum()
    print("\nNormalizing dataset...", end="", flush=True)
    # create a new multiindexed column called rel_freq
    new_col_names = list(processed_barcodes.columns.levels)
    new_col_names[0] = "rel_freq"
    new_cols = pd.MultiIndex.from_product(new_col_names)
    normed_df = processed_barcodes.loc[:, idx["counts"]]/sums
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    # add the new column
    processed_barcodes = pd.concat([processed_barcodes, normed_df], axis=1)

    # calculate medians for rel_freq
    medians = processed_barcodes.loc[idx[:, "WT"], idx["rel_freq"]].median()

    # create a new multiindexed column called rel_wt
    new_col_names = list(processed_barcodes.columns.levels)
    new_col_names[0] = "rel_wt"
    new_cols = pd.MultiIndex.from_product(new_col_names)
    normed_df = processed_barcodes.loc[:, idx["rel_freq"]].divide(medians)
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    # add new column
    processed_barcodes = pd.concat([processed_barcodes, normed_df], axis=1)
    print("Done.", flush=True)
    return processed_barcodes


def regress(df, group_name, experimental_info_csv_path):
    """
    Runs a linear regression across timepoints and returns a dataframe with the slopes
    If average_days is True then the two days will be averaged and one value returned
    instead of two. gen_times is a dictionary mapping days to their  list of generation times
    """
    if group_name in ["APC", "Onion"]:
        warnings.warn("APC and Onion generation times are wonky, aborting...")
        return
    log_df = np.log2(df["rel_wt"])

    try:
        gen_times = utils.load_gen_times(experimental_info_csv_path)[group_name]
    except KeyError as e:
        print("Group not found")
        return

    def _linregress_func(xs, ys):
        slope, intercept, r, p, stderr = scipy.stats.linregress(xs, ys)
        return pd.Series({"slope": slope, "r^2": r**2, "intercept": intercept, "stderr": stderr},
                         index=["slope", "stderr", "r^2", "intercept"], name="stats")

    print("\nRegressing on data...", end="", flush=True)
    df = log_df.stack("days").apply(lambda ys: _linregress_func(gen_times[ys.name[-1]], ys), axis=1)

    print("Done.\n\nThresholding using stderr <0.1 required for all days", flush=True)
    before = len(df.unstack("days"))
    df = df[df["stderr"] < 0.1].unstack("days")
    print("discarding {} barcodes...".format(before - len(df)))
    return df["slope"]


def groupby_filter(df, levels=["codons", "positions"], single_bc_cutoff=0.02, pval_cutoff=0.05):
    """
    Returns a dataframe of groups that pass filtering. Filters using a 2 sample t-test if
    there are more than 1 unique barcodes mapped to a group, otherwise uses
    a distance cutoff. Both filtering cutoffs can be modified by changing the respective
    parameter.

    Parameters:
    df, pd.DataFrame: the data
    levels, list(str): the indices to group together
    """
    if df is None:
        return
    print("\nGrouping and filtering...", end="", flush=True)

    def sig_filter(d1, d2):
        d1, d2 = d1[~pd.isnull(d1)], d2[~pd.isnull(d2)]
        if (len(d1) == 1 or len(d2) == 1) and (d1.mean() - d2.mean())**2 > single_bc_cutoff:
            return pd.Series({"size_d1": len(d1), "size_d2": len(d2)})
        pval = np.nan
        if len(d1) > 1 and len(d2) > 1:
            pval = scipy.stats.ttest_ind(d1, d2)[1]
            if pval < pval_cutoff:
                return pd.Series({"size_d1": len(d1), "size_d2": len(d2), "pval":pval})
        return pd.Series({"mean": pd.concat([d1, d2]).mean(), "size_d1": len(d1), "size_d2": len(d2), "pval": pval})

    df = df.groupby(level=levels).apply(lambda x: sig_filter(x["d1"], x["d2"])).unstack(level=-1)
    print("Done.", flush=True)
    return df.loc[~pd.isnull(df["mean"]), :]

