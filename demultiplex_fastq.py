import re
import pandas as pd
from collections import defaultdict
from sys import argv

_, filename = argv

print(filename)

et0h = {"ATCGTG":"d1t0", "TGAGTG":"d1t1", "CGCCTG":"d1t2", "GCCATG":"d2t0", "AAAATG":"d2t1", "TGTTGG":"d2t2"}
control = {"ATTCCG":"d1t0", "AGCTAG":"d1t1", "GTATAG":"d1t2"}

et0h_barcodes_to_counts = defaultdict(int)
control_barcodes_to_counts = defaultdict(int)

try:
    with open(filename) as lane1_file:
        lines = enumerate(lane1_file)
        for idx, line in lines:
            header = line
            _, seq = next(lines, (None, None))
            next(lines, (None, None))
            _, qual = next(lines, (None, None))
            index = header[-7:].strip()
            barcode = seq[:18]

            if index in et0h:
                et0h_barcodes_to_counts[barcode+"_%s"%et0h[index]] += 1
            elif index in control:
                control_barcodes_to_counts[barcode+"_%s"%control[index]] += 1
except:
    print("Stopping")
finally:
    et0h = pd.Series(et0h_barcodes_to_counts).to_frame().reset_index()
    et0h.columns = ["barcodes_exp", "counts"]
    control = pd.Series(control_barcodes_to_counts).to_frame().reset_index()
    control.columns = ["barcodes_exp", "counts"]

    et0h[["barcodes", "exp"]] = et0h["barcodes_exp"].str.split("_", expand=True)
    et0h.drop("barcodes_exp", axis=1, inplace=True)
    et0h = et0h[["exp", "barcodes", "counts"]]
    et0h.to_csv("et0h_barcodes_to_count.csv", index=False)

    control[["barcodes", "exp"]] = control["barcodes_exp"].str.split("_", expand=True)
    control.drop("barcodes_exp", axis=1, inplace=True)
    control = control[["exp", "barcodes", "counts"]]
    control.to_csv("control_barcodes_to_count.csv", index=False)
