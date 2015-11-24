import pandas as pd
import numpy as np
import sklearn.neighbors
import os
import pickle
from multiprocessing import Pool, Manager, cpu_count
from Bio.Seq import Seq
import time, itertools
from timeit import default_timer as timer


canonical_yeast_ubq = list(" QIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG ")


def load_mutant_map(hdf5_datastorepath, allele_pkl_path = None):

    if allele_pkl_path is None:
        allele_pkl_path = os.path.join(*[os.path.split(hdf5_datastorepath)[0], "allele_dic_with_WT.pkl"])
    with open(allele_pkl_path, "rb") as f:
        barcode_mutant_map = pd.DataFrame(pickle.load(f)).T.reset_index()

    barcode_mutant_map.columns = ["barcodes", "positions", "codons"]
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].astype(np.str)

    barcode_mutant_map["WT"] = barcode_mutant_map["codons"] == "WT"
    # add dummy value for WT barcodes
    barcode_mutant_map["amino acids"] = "WT"
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).transcribe()))
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "amino acids"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).translate()))

    return barcode_mutant_map



def load_gen_times(df, experimental_info_csv_path):
    """
    Loads generation times from the experimental info csv on the Pubs website and
    adds it to the dataframe

    :param df:
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
    primers = primers.set_index(["group", "days", "timepoints"])
    primers = primers.astype(np.float)

    df = df.stack("timepoints").reset_index().set_index(["group", "days", "timepoints"]).join(primers)
    return df.set_index(["barcodes", "codons", "amino acids", "positions"], append=True).unstack("timepoints")

btree = None


def nearest_neighbors(args):
    x, q = args
    dist, ind = btree.query(x, k=2)
    q.put(0)
    return dist, ind


def parallel_hamming(arr, func):
    start = timer()
    num_cpus = cpu_count()
    print("Using {} cpus to compute hamming distances for {} items...\n\nPercent complete: ".format(num_cpus, len(arr)), end="")
    with Pool(num_cpus) as pool:
        m = Manager()
        q = m.Queue()
        result = pool.map_async(func, [(arr[i:i+1000], q) for i in range(0, len(arr), 1000)])
        for i, elem in enumerate(itertools.cycle('\|/-')):
            if result.ready():
                break
            size = q.qsize()
            print("Percent complete: {:.0%} {}".format(size/(len(arr)/1000), elem), end="\r")
            time.sleep(0.25)
        print("Percent complete: {:.0%}".format(1))
    return (np.concatenate(data) for data in zip(*result.get()))


def hamming_correct(raw_barcode_data, mapped_barcode_data, barcode_mutant_map, max_dist=3):
    """
    High performance mapping of unmapped barcodes onto barcode library. Only unambigious and small errors will be corrected.

    :param raw_barcode_data:
    :param mapped_barcode_data:
    :param barcode_mutant_map:
    :param max_dist: The maximum Hamming distance allowed
    :return: new_mapped_barcode_data
    """
    unmapped_raw_barcode_data = raw_barcode_data[~raw_barcode_data["barcodes"].isin(barcode_mutant_map["barcodes"])]
    print("Recovering {} unmapped barcodes".format(len(unmapped_raw_barcode_data)))
    print("Recovering {} counts".format(unmapped_raw_barcode_data["counts"].sum()))

    unmapped_barcodes = unmapped_raw_barcode_data["barcodes"].unique().astype("str")
    unmapped_as_int = unmapped_barcodes.view('S4').reshape((unmapped_barcodes.size, -1)).view(np.uint32)

    barcode_lib = barcode_mutant_map["barcodes"].unique().astype("str")
    barcode_lib_as_int = barcode_lib.view('S4').reshape((barcode_lib.size, -1)).view(np.uint32)

    global btree
    btree = sklearn.neighbors.BallTree(barcode_lib_as_int, leaf_size=40, metric="hamming")

    print("Mapping {} unique unmapped barcodes onto {} library barcodes".format(len(unmapped_barcodes), len(barcode_lib)), flush=True)

    dist, ind = parallel_hamming(unmapped_as_int, nearest_neighbors)

    dist_from_muts = dist*18

    # find only lib barcodes that are nonambiguous and only < max_dist steps away from the
    # unmapped barcode
    mask = (np.diff(dist_from_muts).flatten() >= 1) & (dist_from_muts[:, 0] < max_dist)
    output = ind[mask]
    # preserve index position so we know which unmapped barcode each value
    # corresponds to
    og_idx = np.arange(len(ind))

    corrected = np.vstack([unmapped_barcodes[og_idx[mask]], barcode_lib[output][:, 0]]).T

    print("{} barcodes will be remapped".format(len(corrected)))

    # mapping of the erroneous barcode to its closest neighbor in the barcode lib
    new_mapping = unmapped_raw_barcode_data.merge(pd.DataFrame(corrected, columns=["barcodes", "corrected_bc"]), on="barcodes")
    new_mapping.rename(columns={new_mapping.columns[-1]:"barcodes", new_mapping.columns[-3]:"old_barcodes"}, inplace=True)
    new_mapping.set_index(["group", "days", "timepoints", "barcodes"], inplace=True)
    new_mapping.sort_index(inplace=True)

    # barcodes used to be unique prior to the correction so now group them together because
    # they represent the same barcode
    grouped_new_mapping = new_mapping.groupby(level=[i for i in range(len(new_mapping.index.levels))]).sum()

    mapped_barcode_data_reindexed = mapped_barcode_data.set_index(["group", "days", "timepoints", "barcodes"]).sort_index()
    mapped_barcode_data_reindexed["new_counts"] = mapped_barcode_data_reindexed.loc[:, "counts"].add(grouped_new_mapping["counts"], fill_value=0)
    # len(np.unique(test3.index.get_level_values("barcodes")))
    new_unmapped = grouped_new_mapping[~grouped_new_mapping.index.isin(mapped_barcode_data_reindexed.index)].loc[:, "counts"]

    new_raw_barcode_data = pd.concat([mapped_barcode_data_reindexed["new_counts"], new_unmapped])
    new_raw_barcode_data = new_raw_barcode_data.reset_index(name="counts")

    new_mapped_barcode_data = new_raw_barcode_data.merge(barcode_mutant_map, on="barcodes")

    return new_mapped_barcode_data


def subtract_control(df, merge_on=["amino acids", "positions"]):
    """
    Subtracts control data from the perturbation data

    :param df: either aa_weighted or codons_weighted
    :param merge_on: the columns to merge on either AA-pos or codons-pos
    :return: a dataframe with a difference column
    """
    s1 = df.loc[df.index.get_level_values("group") != "Control", "weighted mean slope"]
    s2 = df.loc["Control", "weighted mean slope"]
    merged = pd.merge(s1.reset_index(), s2.reset_index(), on=merge_on)
    merged["diff"] = merged["weighted mean slope_x"] - merged["weighted mean slope_y"]
    merged.drop("days_y", axis=1, inplace=True)
    merged.rename(columns={'days_x':'days'}, inplace=True)
    merged = merged.set_index(["group", "days", "amino acids", "positions"]).sort_index()
    merged.drop("WT", level="amino acids", inplace=True)
    merged.drop(["weighted mean slope_x", "weighted mean slope_y"], axis=1, inplace=True)
    return merged
