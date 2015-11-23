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


def plot_heatmap(df, title, file_name, vmax=1.5, vmin=-1.5, cmap = sns.diverging_palette(145, 280, s=85, l=25, as_cmap=True)):
    sns.set_style("white")

    fig = plt.figure(figsize=(16, 8))
    df = df.reindex(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'C', 'S', 'T', 'P', 'D', 'E', 'N', 'Q', 'H', 'K', 'R', '*'])
    df.index.values[-1:] = ['STOP']
    df = df.iloc[:, :75]
    sns.heatmap(df, square=True, linewidths=0.25, xticklabels=2, cmap=cmap, vmin=vmin, vmax=vmax)

    plt.yticks(rotation=0)
    plt.title(title, fontsize=16)
    pos = fig.axes[1].get_position()
    fig.axes[1].set_position([pos.x0, 0.4, 0.1, 0.2])
    sns.axlabel("Ub position", "Amino acid")
    plt.savefig(file_name, dpi=300, bbox_inches='tight', pad_inches=0)