import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import utils


def plot_std_vs_log_mean(df, file_name="std_vs_log_mean.png", title="Consistency of mutant fitness increases\nwith average read count", title_fontsize=16, dpi=300):
    """Plots the comparison between log10(read count at t=0) and standard error of per-barcode slopes."""
    out = df.groupby(level=["group", "days", "codons", "positions"]).agg({("fitness", "slope"):np.std, ("counts", "t0"):np.mean})
    out = out[[("counts", "t0"), ("fitness", "slope")]]
    out.columns = ["mean", "std"]
    out["mean"] = np.log(out["mean"])

    sns.set_style("white")
    ax = sns.jointplot("mean", "std", data=out[~pd.isnull(out["std"])], kind="kde", ylim=(0, 0.5), xlim=(0, 10))
    ax.ax_joint.set_xlabel("log-mean of t0 counts per mutant")
    ax.ax_joint.set_ylabel("std of slopes per mutant")
    plt.subplots_adjust(top=0.85)
    plt.suptitle(title, fontsize=title_fontsize, y=0.95)
    if file_name is not None:
        plt.savefig(file_name, dpi=dpi)
    return plt.gcf()


def plot_masked_mutant_wt_heatmap(df, by_codon=False, slope_column='weighted mean slope', filter_fn=lambda row: row['p-value'] > 0.05):
    """A convenience function for calling plot_two_color_masked_heatmap."""
    plot_two_color_masked_heatmap(df, by_codon=by_codon, slope_column=slope_column, filter_fn=filter_fn)


def plot_masked_case_control_heatmap(df, by_codon=False, slope_column='case-control slope diff', filter_fn=lambda row: (row['p-value: case'] > 0.05 and row['p-value: control'] > 0.05) or row['d-value'] > 0.05, weighted=True):
    """A convenience function for calling plot_two_color_masked_heatmap."""
    plot_two_color_masked_heatmap(df, by_codon=by_codon, slope_column=slope_column, filter_fn=filter_fn)


def plot_two_color_masked_heatmap(df, by_codon=False, slope_column='weighted mean slope', filter_fn=None, aggfunc=np.mean,
                                  title=None, title_font_size=30, file_name=None, dpi=300):
    """
    Plots a two-color heatmap of fitnesses indexed by position and AA (or codon), masking cells for which a function is true.
    Especially useful for applying a p-value cutoff.
    slope_column: The name of the value to plot as color (heat)
    filter_fn: Mask every cell for which the function is True
    """
    mutant_type = "codons" if by_codon else "amino acids"
    sns.set_style("white")
    x = df.copy().reset_index()
    if filter_fn is not None:
        x[slope_column] = x.apply(lambda row: float('NaN') if filter_fn(row) else row[slope_column], axis=1) # mask NaNs after
    x = pd.pivot_table(x, index=mutant_type, columns='positions', values=slope_column, aggfunc=aggfunc)
    x = x.reindex(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'C', 'S', 'T', 'P', 'D', 'E', 'N', 'Q', 'H', 'K', 'R', '*'])
    sns.heatmap(x, xticklabels=2, mask=np.isnan(x))
    if title is not None:
        plt.title(title, fontsize=title_font_size)
    if file_name is not None:
        plt.savefig(file_name, dpi=dpi, bbox_inches='tight', pad_inches=0)


def plot_mutant_wt_volcano(x, slope_column='weighted mean slope', log_pval_column='-log10(p-value)', xlim=None, ylim=None,
                      hue=None, size=20, data_point_size=50):
    plot_volcano(x, slope_column=slope_column, log_pval_column=log_pval_column, xlim=xlim, ylim=ylim,
                 hue=hue, size=size, data_point_size=data_point_size)


def plot_case_control_volcano(x, slope_column='case-control slope diff', log_pval_column='-log10(d-value)', xlim=None, ylim=None,
                      hue=None, size=20, data_point_size=50):
    plot_volcano(x, slope_column=slope_column, log_pval_column=log_pval_column, xlim=xlim, ylim=ylim,
                 hue=hue, size=size, data_point_size=data_point_size)


def plot_volcano(x, slope_column, log_pval_column, xlim=None, ylim=None,
                      hue=None, size=20, data_point_size=50):
    """
    Plot a volcano plot.
    Examples: plot_volcano(case_control, slope_column='case-control slope diff', log_pval_column='-log10(d-value)')
              plot_volcano(case_control, slope_column='weighted mean slope', log_pval_column='-log10(p-value)')
    """
    sns.set_style("whitegrid")
    lm = sns.lmplot(data=x, x=slope_column, y=log_pval_column, size=size, fit_reg=False, scatter_kws={"s": data_point_size})
    axes = lm.axes
    if xlim is not None:
        axes[0,0].set_xlim(xlim)
    if ylim is not None:
        axes[0,0].set_ylim(ylim)


def plot_heatmap(df, title, file_name, vmax=1.5, vmin=-1.5, cmap=sns.diverging_palette(220, 20, n=7, as_cmap=True), figsize=(16, 8), dpi=300):

    sns.set_style("white")
    fig = plt.figure(figsize=figsize)
    df = df.reindex(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'C', 'S', 'T', 'P', 'D', 'E', 'N', 'Q', 'H', 'K', 'R', '*'])
    df.index.values[-1:] = ['STOP']
    if 77 in df.columns:
        df = df.loc[:, :76]
    if 0 in df.columns:
        df.drop(0, axis=1, inplace=True)
    sns.heatmap(df, square=True, linewidths=0.25, xticklabels=2, cmap=cmap, vmin=vmin, vmax=vmax)
    yeast_ubq = np.array(utils.canonical_yeast_ubq)
    mask = df.apply(lambda x: x.index.isin(np.where(yeast_ubq == x.name)[0]+1), axis=1)
    df.loc[:, :] = 1
    sns.heatmap(df, mask=~mask, cmap=sns.dark_palette("grey", n_colors=1, as_cmap=True, reverse=True), xticklabels=2, cbar=False, square=True, alpha=0.6)

    plt.yticks(rotation=0)
    plt.title(title, fontsize=16)
    pos = fig.axes[1].get_position()
    fig.axes[1].set_position([pos.x0, 0.4, 0.1, 0.2])
    sns.axlabel("Ub position", "Amino acid")
    if file_name is not None:
        plt.savefig(file_name, dpi=dpi, bbox_inches='tight', pad_inches=0)
    return fig

