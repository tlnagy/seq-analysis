import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_std_vs_log_mean(df):
    out = df.groupby(level=["group", "days", "codons", "positions"]).agg({("fitness", "slope"):np.std, ("counts", "t0"):np.mean})
    out = out[[("counts", "t0"), ("fitness", "slope")]]
    out.columns = ["mean", "std"]
    out["mean"] = np.log(out["mean"])

    sns.set_style("white")
    ax = sns.jointplot("mean", "std", data=out[~pd.isnull(out["std"])], kind="kde", ylim=(0, 0.5), xlim=(0, 10))
    ax.ax_joint.set_xlabel("log-mean of t0 counts per mutant")
    ax.ax_joint.set_ylabel("std of slopes per mutant")
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Consistency of mutant fitness increases\nwith average read count", fontsize=16, y=0.95)
    plt.savefig("std_vs_log_mean.png", dpi=300)
