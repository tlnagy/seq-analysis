import pandas as pd
import numpy as np


def load_gen_times(experimental_info_csv_path):
    """
    Loads generation times from the experimental info csv on the Pubs website

    :param experimental_info_csv_path:
    :return:
    """
    primers = pd.read_csv(experimental_info_csv_path)
    primers["exp"] = primers["Sample"].str.replace("Day", "d").str.replace(" T=", "t").str.split(" ").str.get(1)
    primers["group"] = primers["Sample"].str.split(" ").str.get(0)
    primers = primers[["group", "exp", "Generations"]]
    primers = primers[~pd.isnull(primers["group"])].fillna(0).replace("-", 0)
    primers[["days", "timepoints"]] = primers["exp"].str.replace(r"t", "_t").str.split("_", expand=True)
    primers = primers.drop("exp", axis=1)
    gen_times = {}
    for group, group_data in primers.groupby("group"):
        gen_times.update({group: {day: data["Generations"].astype(np.float).tolist() for day, data in group_data.groupby("days")}})
    return gen_times
