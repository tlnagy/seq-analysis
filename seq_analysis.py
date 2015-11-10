from __future__ import print_function, division
import warnings
import sys
try:
    import platform
    if platform.python_version().startswith('2'):
        warnings.warn("Python 2 is old and yucky, please use 3", RuntimeWarning)
except ImportError:
    sys.exit(2)
from Bio.Seq import Seq
import pandas as pd
import pickle
import numpy as np
import scipy.stats


def map_dataset(csv_filepath):
    """
    Loads raw demultiplexed read csv and maps it onto the barcode-mutant
    pairing. Does not process the data further.
    """
    raw_barcode_data = pd.read_csv(csv_filepath)
    num_unique_reads = len(raw_barcode_data)
    unique_read_count_total = raw_barcode_data["counts"].sum()

    with open("allele_dic_with_WT.pkl", "rb") as f:
        barcode_mutant_map = pd.DataFrame(pickle.load(f)).T.reset_index()

    barcode_mutant_map.columns = ["barcodes", "positions", "codons"]
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].astype(np.str)

    # split out days and timepoints into separate columns
    raw_barcode_data[["days", "timepoints"]] = raw_barcode_data["exp"].str.replace(r"t", "_t").str.split("_", expand=True)

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
    for name, group in mapped_barcode_data.groupby("exp"):
            print("{}    {:.2%}".format(name, len(group)/len(barcode_mutant_map)))

    return raw_barcode_data, mapped_barcode_data


def process_data(mapped_barcode_data):
    # transform data
    processed_barcodes = mapped_barcode_data.pivot_table(index=["days", "timepoints", "barcodes", "codons","amino acids", "positions"],
                                                         values=["counts", "rel_freq", "rel_wt"])

    processed_barcodes = processed_barcodes.unstack("days").unstack("timepoints")

    idx = pd.IndexSlice
    # Throw out barcodes that have no counts in at least one experiment
    before = len(processed_barcodes)
    processed_barcodes = processed_barcodes[pd.notnull(processed_barcodes.loc[:, idx["counts"]]).sum(axis=1) == len(processed_barcodes.columns)]
    print("\nDiscarding {} barcodes that had a count of 0 in at least one timepoint".format(before-len(processed_barcodes)))

    before = len(processed_barcodes)
    # Throw out values that have a count 1 in any experiment
    processed_barcodes = processed_barcodes[processed_barcodes[processed_barcodes == 1].count(axis=1) == 0]
    print("\nDiscarding {} barcodes that had a count of 1 at any timepoint".format(before-len(processed_barcodes)))

    sums = processed_barcodes["counts"].sum()

    # create a new multiindexed column called rel_freq
    new_col_names = list(processed_barcodes.columns.levels)
    new_col_names[0] = "rel_freq"
    new_cols = pd.MultiIndex.from_product(new_col_names)
    normed_df = processed_barcodes.loc[:, idx["counts"]]/sums
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    # add the new column
    processed_barcodes = pd.concat([processed_barcodes, normed_df], axis=1)

    #calculate medians for rel_freq
    medians = processed_barcodes.loc[idx[:, "WT"], idx["rel_freq"]].median()

    # create a new multiindexed column called rel_wt
    new_col_names = list(processed_barcodes.columns.levels)
    new_col_names[0] = "rel_wt"
    new_cols = pd.MultiIndex.from_product(new_col_names)
    normed_df = processed_barcodes.loc[:, idx["rel_freq"]].divide(medians)
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    #add new column
    processed_barcodes = pd.concat([processed_barcodes, normed_df], axis=1)
    return processed_barcodes


def _linregress_func(xs, ys):
    slope, intercept, r, p, stderr = scipy.stats.linregress(xs, ys)
    return pd.Series({"slope": slope, "r^2": r**2, "intercept": intercept, "stderr": stderr},
                     index=["slope", "stderr", "r^2", "intercept"], name="stats")


def regress(df, gen_times={"d1": [0, 3.14, 5.14], "d2": [0, 1.76, 4.02]}):
    """
    Runs a linear regression across timepoints and returns a dataframe with the slopes
    If average_days is True then the two days will be averaged and one value returned
    instead of two. gen_times is a dictionary mapping days to their  list of generation times
    """
    log_df = np.log2(df["rel_wt"])

    print("regressing...\n", flush=True)
    df = log_df.stack("days").apply(lambda ys: _linregress_func(gen_times[ys.name[-1]], ys), axis=1)

    print("thresholding using stderr <0.1 required for all days")
    df = df.unstack("days")
    before = len(df)
    df = df[df[df["stderr"] < 0.1]["stderr"].count(axis=1) == len(df["stderr"].columns)]
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

    def sig_filter(d1, d2):
        if len(d1) == 1 and (float(d1) - float(d2))**2 > single_bc_cutoff:
            return pd.Series({"size":len(d1)})
        pval = np.nan
        if len(d1) > 1:
            pval = scipy.stats.ttest_ind(d1, d2)[1]
            if pval < pval_cutoff:
                return pd.Series({"size":len(d1), "pval":pval})
        return pd.Series({"mean":pd.concat([d1, d2]).mean(), "size":len(d1), "pval":pval})

    df = df.groupby(level=levels).apply(lambda x: sig_filter(x["d1"], x["d2"])).unstack(level=-1)
    return df.loc[~pd.isnull(df["mean"]), :]

