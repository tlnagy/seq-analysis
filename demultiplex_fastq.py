from collections import defaultdict
import pandas as pd
import numpy as np
from sys import argv

_, filename = argv

primers = pd.read_csv("truseq_primers.csv")
primers["exp"] = primers["Sample"].str.replace("Day", "d").str.replace(" T=", "t").str.split(" ").str.get(1)
primers["group"] = primers["Sample"].str.split(" ").str.get(0)
primers = primers[["group", "exp", "Indexing BC"]]

group_exp_barcode_map = primers[~pd.isnull(primers["group"])].set_index("Indexing BC").to_dict()
exp, group = group_exp_barcode_map["exp"], group_exp_barcode_map["group"]
data = {key:defaultdict(int) for key in group}

with open(filename) as fastq_file, open(filename+".log", "w") as logfile:
    lines = enumerate(fastq_file)
    for idx, line in lines:
        header = line
        _, seq = next(lines, (None, None))
        next(lines, (None, None))
        _, qual = next(lines, (None, None))
        index = header[-7:].strip()
        barcode = seq[:18]

        if index in group:
            print(group[index], file=logfile)
            data[index][barcode] += 1

df = pd.DataFrame([exp, group, data]).T
df.columns = ["exp", "group", "data"]
df2 = df["data"].apply(lambda x: pd.Series(x))
df2.columns.name = "barcodes"
df2.index.name = "index"
df.drop("data", axis=1, inplace=True)
results = df.reset_index().merge(df2.stack("barcodes").reset_index(), on="index").set_index(["group", "index", "exp", "barcodes"])
results.columns = ["counts"]
results.reset_index(level="exp", inplace=True)
results[["days", "timepoints"]] = results["exp"].str.replace(r"t", "_t").str.split("_", expand=True)
results = results.reset_index().drop(["index", "exp"], axis=1)
str_cols = ["group", "days", "timepoints", "barcodes"]
results[str_cols] = results[str_cols].astype(np.str)
results.set_index(str_cols, inplace=True)
print(results)
results.to_hdf(filename+".h5", key="raw_barcode_data", mode="w", complevel=9)

