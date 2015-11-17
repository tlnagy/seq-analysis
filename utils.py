import pandas as pd
import numpy as np


def load_gen_times(df, experimental_info_csv_path):
    """
    Loads generation times from the experimental info csv on the Pubs website and
    adds it to the dataframe

    :param df:
    :param experimental_info_csv_path:
    :return:
    """
    primers = pd.read_csv("seq-analysis/truseq_primers.csv")
    primers["exp"] = primers["Sample"].str.replace("Day", "d").str.replace(" T=", "t").str.split(" ").str.get(1)
    primers["group"] = primers["Sample"].str.split(" ").str.get(0)
    primers = primers[["group", "exp", "Generations"]]
    primers = primers[~pd.isnull(primers["group"])].fillna(0).replace("-", 0)
    primers[["days", "timepoints"]] = primers["exp"].str.replace(r"t", "_t").str.split("_", expand=True)
    primers = primers.drop("exp", axis=1)
    primers = primers.set_index(["group", "days", "timepoints"])
    primers = primers.astype(np.float)

    df = df.stack("timepoints").reset_index().set_index(["group", "days", "timepoints"]).join(primers)
    return df.set_index(["barcodes", "codons", "amino acids", "positions"], append=True).unstack("timepoints")
