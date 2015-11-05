from Bio.Seq import Seq
import pandas as pd
import pickle
import numpy as np

def load_dataset(csv_filepath):
    barcode_counts = pd.read_csv(csv_filepath)
    num_unique_reads = len(barcode_counts)
    unique_read_count_total = barcode_counts["counts"].sum()

    data = pickle.load(open("allele_dic_with_WT.pkl", "rb"))

    data = pd.DataFrame(data).T.reset_index()
    data.columns = ["barcodes", "positions", "codons"]
    data["barcodes"] = data["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    data["barcodes"] = data["barcodes"].astype(np.str)
    barcode_counts[["days","timepoints"]] = barcode_counts["exp"].str.replace(r"t", "_t").str.split("_", expand=True)

    data["WT"] = data["codons"] == "WT"
    data["amino acids"] = "@"
    data.loc[~data["WT"], "codons"] = data.loc[~data["WT"], "codons"].apply(lambda x: str(Seq(x).transcribe()))
    data.loc[~data["WT"], "amino acids"] = data.loc[~data["WT"], "codons"].apply(lambda x: str(Seq(x).translate()))

    barcode_counts = barcode_counts.merge(data, on="barcodes", how="inner")
    barcode_counts["positions"] = barcode_counts["positions"].astype(np.int)
    barcode_counts.index.name = "idxs"

    print("Percent of unique barcodes mapped: {:.2%}".format(len(barcode_counts)/num_unique_reads))
    print("Percent of total reads mapped: {:.2%}\n".format(barcode_counts["counts"].sum()/unique_read_count_total))

    print("group   % of total barcodes found")
    for name, group in barcode_counts.groupby("exp"):
        print("{}    {:.2%}".format(name, len(group)/len(data)))

    # transform data
    df = barcode_counts.pivot_table(index=["days", "timepoints", "barcodes", "codons","amino acids", "positions"],  
                                         values=["counts"])

    df = df.unstack("days").unstack("timepoints")

    idx = pd.IndexSlice
    # Throw out values that have no counts in at least one 
    df = df[pd.notnull(df.loc[:, idx["counts"]]).sum(axis=1) == 6]
    sums = df["counts"].sum()
    df = df.stack("days").stack("timepoints")
    flattened_df = df.unstack("days").unstack("timepoints")
    df.loc[:, "rel_freq"] = (flattened_df.loc[:,
        idx["counts"]]/sums).stack("days").stack("timepoints")
    medians = df.loc[idx[:, "WT"], idx["rel_freq"]].median()
    flattened_df = df.unstack("days").unstack("timepoints")
    df.loc[:, "rel_wt"] = (flattened_df.loc[:,
        idx["rel_freq"]]/medians).stack("days").stack("timepoints")
    df = df.unstack("days").unstack("timepoints")
    return(df)
