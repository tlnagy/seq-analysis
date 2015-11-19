from __future__ import print_function, division
import warnings
import sys, os
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
import utils
from multiprocessing import Pool, Manager, cpu_count
import itertools, time
from timeit import default_timer as timer


def process_data(hdf5_datastorepath, allele_pkl_path = None, experimental_info_csv_path=None, hamming_correct=False):
    print("Loading data...", flush=True, end="")
    idx = pd.IndexSlice
    raw_barcode_data = pd.read_hdf(hdf5_datastorepath, key="grouped_data")

    if allele_pkl_path is None:
        allele_pkl_path = os.path.join(*[os.path.split(hdf5_datastorepath)[0], "allele_dic_with_WT.pkl"])
    with open(allele_pkl_path, "rb") as f:
        barcode_mutant_map = pd.DataFrame(pickle.load(f)).T.reset_index()
    print("Done.\n\nMapping...", flush=True, end="")

    barcode_mutant_map.columns = ["barcodes", "positions", "codons"]
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].apply(lambda x: str(Seq(x).reverse_complement()))
    barcode_mutant_map["barcodes"] = barcode_mutant_map["barcodes"].astype(np.str)
    raw_barcode_data = raw_barcode_data.reset_index()

    barcode_mutant_map["WT"] = barcode_mutant_map["codons"] == "WT"
    # add dummy value for WT barcodes
    barcode_mutant_map["amino acids"] = "WT"
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).transcribe()))
    barcode_mutant_map.loc[~barcode_mutant_map["WT"], "amino acids"] = barcode_mutant_map.loc[~barcode_mutant_map["WT"], "codons"].apply(lambda x: str(Seq(x).translate()))

    mapped_barcode_data = raw_barcode_data.merge(barcode_mutant_map, on="barcodes", copy=False)
    print("Done.", flush=True)

    if hamming_correct:
        mapped_barcode_data = utils.hamming_correct(raw_barcode_data, mapped_barcode_data, barcode_mutant_map)

    processed_barcodes = mapped_barcode_data.pivot_table(index=["group", "days", "timepoints", "barcodes", "codons","amino acids", "positions"],
                                                         values=["counts"])
    processed_barcodes = processed_barcodes.unstack("timepoints")

    # Onion d1t2 is simply a repeat of t1 so we'll toss it
    processed_barcodes.loc[idx["Onion", "d1"], idx["counts", "t2"]] = np.nan

    # mapping transitions between timepoints, -1 indicates True to False aka NaN to not-NaN which should be eliminated
    # checks if any such transition occurs in a row
    has_nan_to_val_trans = np.any(np.diff(pd.isnull(processed_barcodes).values.view(np.int8)) == -1, axis=1)
    print("Throwing out barcodes that have nans prior to counts...\n{}".format(processed_barcodes[has_nan_to_val_trans].groupby(level="group").size()))
    processed_barcodes = processed_barcodes[~has_nan_to_val_trans]

    sums = processed_barcodes.unstack("group").unstack("days").reorder_levels([0, 2, 3, 1], axis=1).sort_index(level=["group", "days"], axis=1).sum()

    rel_freq = processed_barcodes.unstack("group").unstack("days").reorder_levels([0, 2, 3, 1], axis=1).sort_index(level=["group", "days"], axis=1)/sums
    rel_freq = rel_freq.stack("group").stack("days").reorder_levels([4, 5, 0, 1, 2, 3]).sort_index(level=["group", "days"])

    medians = rel_freq.unstack("group").unstack("days").loc[idx[:, "WT"], :].median()

    rel_wt = rel_freq.unstack("group").unstack("days")/medians

    rel_wt = rel_wt.stack("group").stack("days").reorder_levels([4, 5, 0, 1, 2, 3]).sort_index(level=["group", "days"])
    rel_wt.columns.set_levels(["rel_wt"], level=0, inplace=True)
    rel_wt = np.log2(rel_wt)

    if experimental_info_csv_path is None:
        experimental_info_csv_path = os.path.join(*[os.path.split(hdf5_datastorepath)[0], "truseq_primers.csv"])
    rel_wt_gen_times = utils.load_gen_times(rel_wt, experimental_info_csv_path)

    all_tp = rel_wt_gen_times[np.all(pd.notnull(rel_wt_gen_times["rel_wt"]), axis=1)]

    slopes = groupby_parallel(all_tp.groupby(level=["group", "days", "amino acids"]), linregress_df)

    # handle barcodes that only show up in two timepoints
    two_tp_only = rel_wt_gen_times[pd.isnull(rel_wt_gen_times["rel_wt"]).sum(axis=1) == 1]
    deltas = two_tp_only.diff(axis=1)
    two_tp_slopes = (deltas[("rel_wt", "t1")]/deltas[("Generations", "t1")])
    two_tp_slopes.name = "slope"
    slopes = pd.concat([slopes, two_tp_slopes.to_frame()])
    slopes.sort_index(inplace=True)

    # add counts data back onto the slope dataframe
    slopes["data"] = "fitness"
    slopes.set_index("data", append=True, inplace=True)
    slopes = slopes.unstack("data").reorder_levels([1, 0], axis=1)
    new_cols = pd.MultiIndex.from_product(processed_barcodes.columns.levels)
    normed_df = processed_barcodes.loc[processed_barcodes.index.isin(slopes.index), idx["counts"]]
    normed_df = pd.DataFrame(normed_df.values, index=normed_df.index, columns=new_cols)
    slopes = pd.concat([slopes, normed_df], axis=1)

    extra_funcs = {("fitness", "slope"):[len, np.std], ("counts", "t0"):np.sum}
    weighted_avg_func = lambda x: (x[("fitness", "slope")]*np.log(x[("counts", "t0")])/np.log(x[("counts", "t0")]).sum()).sum()
    extra_col_names = ["# unique barcodes", "stddev of slope", "sum of t0 reads"]

    aa_weighted = slopes.groupby(level=["group", "days", "amino acids", "positions"]).apply(weighted_avg_func)
    aa_weighted.name = "weighted mean slope"
    aa_extra = slopes.groupby(level=["group", "days", "amino acids", "positions"]).agg(extra_funcs)
    aa_extra.columns = extra_col_names
    aa_weighted = pd.concat([aa_weighted, aa_extra], axis=1)
    codons_weighted = slopes.groupby(level=["group", "days", "codons", "positions"]).apply(weighted_avg_func)
    codons_weighted.name = "weighted mean slope"
    codons_extra = slopes.groupby(level=["group", "days", "codons", "positions"]).agg(extra_funcs)
    codons_extra.columns = extra_col_names
    codons_weighted = pd.concat([codons_weighted, codons_extra], axis=1)

    return slopes, codons_weighted, aa_weighted


def linregress_wrapper(xs, ys):
    slope, intercept, r, p, stderr = scipy.stats.linregress(xs, ys)
    return pd.Series({"slope": slope, "r^2": r**2, "intercept": intercept, "stderr": stderr},
                     index=["slope", "intercept", "stderr"], name="stats")


def linregress_df(args):
    df, q = args
    result = df.apply(lambda x: linregress_wrapper(x.values[3:], x.values[:3]), axis=1)
    q.put(0)
    return result


def groupby_parallel(groupby_df, func):
    start = timer()
    num_cpus = cpu_count()
    print("\nUsing {} cpus in parallel...\n\nPercent complete: ".format(num_cpus), end="")
    with Pool(num_cpus) as pool:
        m = Manager()
        q = m.Queue()
        result = pool.map_async(func, [(group, q) for name, group in groupby_df])
        for i, elem in enumerate(itertools.cycle('\|/-')):
            if result.ready():
                break
            size = q.qsize()
            print("Percent complete: {:.0%} {}".format(size/len(groupby_df), elem), end="\r")
            time.sleep(0.25)
        print("Percent complete: {:.0%}".format(1))
    print("\nProcessed on {} slopes in {:.1f}s".format(groupby_df.count().sum().max(), timer() - start))
    return pd.concat(result.get())


