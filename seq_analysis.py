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

raw_barcode_data = None
mapped_barcode_data = None


def map_dataset(csv_filepath):
    """
    Loads raw demultiplexed read csv and maps it onto the barcode-mutant
    pairing. Does not process the data further.
    """
    global raw_barcode_data
    raw_barcode_data = pd.read_csv(csv_filepath)
    num_unique_reads = len(raw_barcode_data)
    unique_read_count_total = raw_barcode_data["counts"].sum()

    with open("allele_dic_with_WT.pkl", "rb") as f:
        barcode_mutant_map = pd.DataFrame(pickle.load(f)).T.reset_index()

    barcode_mutant_map.columns = ["barcodes", "positions", "codons"]
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].astype(np.str)

    # split out days and timepoints into separate columns
    raw_barcode_data[["days","timepoints"]] = raw_barcode_data["exp"].str.replace(r"t", "_t").str.split("_", expand=True)

    barcode_mutant_map["WT"] = barcode_mutant_map["codons"] == "WT"
    # add dummy value for WT barcodes
    barcode_mutant_map["amino acids"] = "@"
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).transcribe()))
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "amino acids"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).translate()))

    global mapped_barcode_data
    mapped_barcode_data = raw_barcode_data.merge(barcode_mutant_map, on="barcodes", how="inner")
    mapped_barcode_data["positions"] = mapped_barcode_data["positions"].astype(np.int)
    mapped_barcode_data.index.name = "idxs"

    print("Percent of unique barcodes mapped: {:.2%}".format(len(mapped_barcode_data)/num_unique_reads))
    print("Percent of total reads mapped: {:.2%}\n".format(mapped_barcode_data["counts"].sum()/unique_read_count_total))

    print("group   % of total barcodes found")
    for name, group in mapped_barcode_data.groupby("exp"):
            print("{}    {:.2%}".format(name, len(group)/len(barcode_mutant_map)))

    return raw_barcode_data, mapped_barcode_data


def process_data(csv_filepath = "et0h_barcodes_to_count.csv"):
    if mapped_barcode_data is None:
        print("Data not yet loaded, loading now...\n")
        map_dataset(csv_filepath)

    # transform data
    df = mapped_barcode_data.pivot_table(index=["days", "timepoints", "barcodes", "codons","amino acids", "positions"],
                                         values=["counts", "rel_freq", "rel_wt"])

    df = df.unstack("days").unstack("timepoints")

    idx = pd.IndexSlice
    # Throw out barcodes that have no counts in at least one experiment
    df = df[pd.notnull(df.loc[:, idx["counts"]]).sum(axis=1) == 6]

    before = len(df)
    # Throw out values that have a count 1 in any experiment
    df = df[df[df == 1].count(axis=1) == 0]
    print("\nDiscarding {} barcodes that had a count of 1 at any timepoint".format(before-len(df)))

    sums = df["counts"].sum()

    # create a new multiindexed column called rel_freq
    new_cols = pd.MultiIndex.from_product([["rel_freq"], ["d1", "d2"], ["t0", "t1", "t2"]])
    normed_df = df.loc[:, idx["counts"]]/sums
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    # add the new column
    df = pd.concat([df, normed_df], axis=1)

    #calculate medians for rel_freq
    medians = df.loc[idx[:, "WT"], idx["rel_freq"]].median()
    median_string = medians.to_string(float_format='{:.8g}'.format)
    print("\nMedian values for relative frequencies\n{}".format(median_string))

    # create a new multiindexed column called rel_wt
    normed_df = df.loc[:, idx["rel_freq"]].divide(medians)
    new_cols = pd.MultiIndex.from_product([["rel_wt"], ["d1", "d2"], ["t0", "t1", "t2"]])
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)

    #add new column
    df = pd.concat([df, normed_df], axis=1)
    return df
