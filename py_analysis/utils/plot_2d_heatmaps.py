import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_probability_heatmaps(
    df,
    id_col: str = 'id',
    x_col: str = 'left_open',
    y_col: str = 'right_open',
    value_col: str = 'df_exp_norm',
    log_scale: bool = False,
    colormap: str = 'plasma_r',
    vmin: float = None,
    vmax: float = None,
    figsize: tuple = None,      # e.g. (4*n_panels, 6)
    sharex: bool = True,
    sharey: bool = False,
    constrained_layout: bool = True,
    dpi: int = 150,
    xlabel: str = None,
    ylabel: str = None,
    colorbar_label: str = None,
    suptitle: str = None,
    title_fontsize: float = 14,
    panel_title_fontsize: float = 12,
    label_fontsize: float = 10,
    tick_fontsize: float = 8,
    xtick_rotation: float = 45
):
    """
    Plot one panel per unique id showing a heatmap of `value_col` over (x_col, y_col).

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns [id_col, x_col, y_col, value_col].
    id_col : str
        Column to split panels by (slider/grouping).
    x_col, y_col : str
        Names of the DataFrame columns for the horizontal (x) and vertical (y) axes.
    value_col : str
        Name of the column to color by.
    log_scale : bool
        If True, use a logarithmic color scale (LogNorm). Otherwise linear.
    colormap : str
        Any Matplotlib colormap name.
    vmin, vmax : float or None
        Color‐scale bounds. If None, computed from the data (nonzero‐min for log).
    figsize : (w, h) or None
        Figure size. Defaults to (4*n, 6).
    sharex, sharey : bool
        Whether to share x/y axes across panels.
    constrained_layout : bool
        Whether to use constrained_layout for nicer spacing.
    dpi : int
        Figure DPI.
    xlabel, ylabel : str or None
        Axis titles (default to x_col, y_col).
    colorbar_label : str or None
        If None, defaults to `value_col + (' (log)' or '')`.
    suptitle : str or None
        Overall title. If None, defaults to `"Probability Heatmaps Across IDs"`.
    title_fontsize, panel_title_fontsize, label_fontsize, tick_fontsize, xtick_rotation
        Font sizes & tick rotation for panels.

    Returns
    -------
    fig, axs : matplotlib Figure and Axes array
    """
    # prepare panel IDs
    unique_ids = np.sort(df[id_col].unique())
    n = len(unique_ids)
    # defaults
    if figsize is None:
        figsize = (4*n, 6)
    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col
    if colorbar_label is None:
        suffix = ' (log scale)' if log_scale else ''
        colorbar_label = f'{value_col}{suffix}'
    if suptitle is None:
        suptitle = "Probability Heatmaps Across IDs"
    # color‐scale bounds
    data_vals = df[value_col]
    if log_scale:
        nz = data_vals[data_vals > 0]
        vmin_ = nz.min() if vmin is None else vmin
        vmax_ = data_vals.max()   if vmax is None else vmax
        norm = LogNorm(vmin=vmin_, vmax=vmax_)
    else:
        vmin_ = data_vals.min() if vmin is None else vmin
        vmax_ = data_vals.max() if vmax is None else vmax
        norm = None

    # create subplots
    fig, axs = plt.subplots(
        1, n, figsize=figsize,
        sharex=sharex, sharey=sharey,
        constrained_layout=constrained_layout,
        dpi=dpi
    )
    if n == 1:
        axs = [axs]

    # plot each panel
    for ax, id_val in zip(axs, unique_ids):
        df_id = df[df[id_col] == id_val]
        pivot = df_id.pivot(
            index=y_col,
            columns=x_col,
            values=value_col
        )
        
        if log_scale:
            # only pass norm, no vmin/vmax
            im = ax.imshow(
                pivot.values,
                origin='lower',
                aspect='auto',
                cmap=colormap,
                norm=LogNorm(vmin=vmin_, vmax=vmax_)
            )
        else:
            # only pass vmin/vmax, no norm
            im = ax.imshow(
                pivot.values,
                origin='lower',
                aspect='auto',
                cmap=colormap,
                vmin=vmin_,
                vmax=vmax_
            )


        # ticks/labels only bottom & left
        ax.tick_params(axis='x', top=False, labeltop=False)
        ax.tick_params(axis='y', right=False, labelright=False)

        ax.set_title(f"ID {id_val}", fontsize=panel_title_fontsize, pad=8)
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.set_xticks(np.arange(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=xtick_rotation,
                           ha='right', fontsize=tick_fontsize)
        ax.set_yticks(np.arange(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=tick_fontsize)

    # shared colorbar
    cbar = fig.colorbar(
        im, ax=axs, orientation='vertical',
        fraction=0.02, pad=0.04
    )
    cbar.set_label(colorbar_label, fontsize=label_fontsize)

    # overall title
    fig.suptitle(suptitle, fontsize=title_fontsize, y=1.02)

    return fig, axs
