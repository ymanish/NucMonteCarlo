from typing import Tuple, List, Sequence
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from py_analysis.config.eads_var import PHOSPHATE_SITES
from py_analysis.config.custom_types import FitResult, EadsResult, Orientation
from py_analysis.adsorbpipe.helpers.fit_multigaussian import _multi_gaussian

##############PLOTING FUNCTIONS FOR RAW DWELL TIME######################

def plot_rawdwell_times(df_left: pd.DataFrame, 
                    df_right: pd.DataFrame,
                    phosphate_sites: List[int] = PHOSPHATE_SITES,
                    title: str = "Dwelltime vs. bp_pos",
                    figsize: Tuple[int, int] = (10, 6)) -> None:
    """Plot dwell times with phosphate site markers."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(df_left['bp_pos'], df_left['dwelltime'], lw=2, color='black', label='Left Zip')
    ax.plot(df_right['bp_pos'], df_right['dwelltime'], lw=2, color='red', label='Right Zip')
    
    for site in phosphate_sites:
        ax.axvline(site, color='gray', linestyle='--', lw=1)
    
    ax.set_xlabel('bp_pos')
    ax.set_ylabel('dwelltime')
    ax.set_title(title)
    ax.legend()
    plt.show()
    # fig.tight_layout()
    # plt.pause(show_time_sec)
    # plt.close(fig)
    return None

##############PLOTING FUNCTIONS FOR CUMULATIVE DWELL TIME######################


def plot_cumulative_dwell_times(left_cum: pd.DataFrame, right_cum: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(left_cum['bp_pos'], left_cum['cumulative_dwell_time'], lw=2, color='blue', label='Left Cumulative')
    ax.plot(right_cum['bp_pos'], right_cum['cumulative_dwell_time'], lw=2, color='orange', label='Right Cumulative')
    for site in PHOSPHATE_SITES:
        ax.axvline(site, color='gray', linestyle='--', lw=1)
    ax.set(xlabel='bp_pos', ylabel='Cumulative Dwell Time', title='Cumulative Dwell Time vs. bp_pos')
    ax.legend()
    plt.show()
    return None





#############PLOTING FUNCTIONS FOR GAUSSIAN FITTING######################


def plot_detected_peaks(dwelltime: np.ndarray | pd.Series,
                        bp_pos: np.ndarray,
                        peaks: Sequence[int],
                        *,
                        title: str, ax: plt.Axes = None) -> None:
    """Quick look at the raw dwell time trace + detected peaks."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(bp_pos, dwelltime, "k-", lw=1.5, label="trace")
    ax.scatter(peaks, dwelltime[np.isin(bp_pos, peaks)], c="red", zorder=3, label="peaks")
    for site in PHOSPHATE_SITES:
        ax.axvline(site, color="grey", lw=1, ls="--")
    ax.set_xlabel("bp position")
    ax.set_ylabel("dwell time (s)")
    ax.set_title(title)
    ax.legend()
    ax.figure.tight_layout()
    return None


def plot_gaussian_fit(dwelltime: np.ndarray | pd.Series,
                        bp: np.ndarray,
                        fit_result: FitResult,
                        *,
                        title: str,
                        show_individual: bool = True, ax: plt.Axes = None) -> None:
    
    """Plot data, composite fit, and optionally individual Gaussians."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(bp, dwelltime, "k", lw=1.5, label="data")
    ax.plot(bp, _multi_gaussian(bp, *fit_result.pars), lw=2, label="composite")

    if show_individual:
        for A, mu, sigma in zip(fit_result.amplitudes, fit_result.mus, fit_result.sigmas):
            g = A * np.exp(-(bp - mu) ** 2 / (2.0 * sigma ** 2))
            ax.plot(bp, g, ls="--", lw=1.2)

    ax.vlines(fit_result.centres, 0, max(dwelltime) * 1.05, color="lightgrey", ls=":", lw=0.8)
    ax.set_title(title)
    ax.set_xlabel("bp")
    ax.set_ylabel("dwell time (s)")
    ax.legend()
    ax.figure.tight_layout()
    return None


def plot_shifted_peaks(dwelltime: np.ndarray | pd.Series,
                        bp: np.ndarray,
                        detected_peaks: Sequence[int],
                        shifted_sites: Sequence[int],
                        values: Sequence[float],  ### These values are same as the original dwell times, using just a reference to show the shifted peaks
                        *,
                        offset: float = 0.0,
                        title: str, ax: plt.Axes = None) -> None:
    """Red dots = detected peaks; blue stars = where they ended up."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(bp, dwelltime, "k", lw=1.5) ### This is the original dwell time trace
    ax.scatter(np.array(detected_peaks), dwelltime[np.isin(bp, detected_peaks)],
                c="red", zorder=3, s=50, label="detected")
    
    ax.scatter(shifted_sites, values,
                marker="*", c="blue", zorder=4,  s=50, label="assigned")
    

    #  Offset arrow (annotation only)
    if offset > 0:
        y_arrow = dwelltime.max() * 0.95
        start_x = bp.min() + 5
        length  = 15  
        end_x   = start_x + length
        ax.annotate("", xy=(end_x, y_arrow), xytext=(start_x, y_arrow),
            arrowprops=dict(arrowstyle="->", lw=2, color="red"))

    elif offset < 0:
        y_arrow = dwelltime.max() * 0.95
        start_x = bp.max() - 5
        length  = 15
        end_x   = start_x -length
        ax.annotate("", xy=(end_x, y_arrow), xytext=(start_x, y_arrow),
                    arrowprops=dict(arrowstyle="->", lw=2, color="red"))


    for site in PHOSPHATE_SITES:
        ax.axvline(site, color="grey", ls="--", lw=1)
    ax.set_xlabel("bp position")
    ax.set_ylabel("dwell / amplitude (s)")
    ax.set_title(title)
    ax.legend()
    ax.figure.tight_layout()

    return None


######### PLOTING FUNCTIONS FOR EADS AND DELTA_G ######################


def plot_tau_dg(results:list[EadsResult], seq_id:str, *, figsize=(12, 8)):
    """
    Plot dwell time (tau) and step ΔG for left/right unwrapping
    in **original** and **flipped** orientations.

    Parameters
    ----------
    results : list[EadsResult]
        Typically '[eads_original, eads_flipped]'.  Order decides
        which row is plotted first.
    figsize : tuple[int, int]
        Size passed to 'plt.subplots()'.

    Notes
    -----
    • Solid line with circles  = τ  (left y-axis)  
    • Dashed line with crosses = ΔG (right y-axis)
    • Plots with same color are same sequences 
    """
    if len(results) != 2:
        raise ValueError("Provide exactly two EadsResult objects (orig and flip)")

    # ensure deterministic row order: ORIGINAL row 0, FLIPPED row 1
    res_sorted = sorted(results, key=lambda r: r.orientation.value)

    fig, axs = plt.subplots(2, 2, figsize=figsize, sharex=True)
    fig.tight_layout(pad=2.0)

    # colour matrix -> diagonal match
    c_diag, c_off = "tab:blue", "tab:orange"
    colour = [[c_diag, c_off],        
              [c_off, c_diag]]      

    # --- helper for one axis --------------------------------------------------
    def _axis(ax, x, tau, dG, side_label, color):
        ax.plot(x, tau, marker="o", linestyle="-", color=color)
        ax.set_ylabel(f"τ {side_label}", color=color)
        ax.tick_params(axis="y", labelcolor=color)
        ax.grid(True, linestyle="-", alpha=0.6)
        sns.despine(ax=ax, left=True, bottom=True)

        ax2 = ax.twinx()
        ax2.plot(x, dG, marker="x", linestyle="--", color=color)
        ax2.set_ylabel(f"ΔG {side_label}", color=color)
        ax2.tick_params(axis="y", labelcolor=color)
        ax2.grid(True, linestyle="--", alpha=0.6)


    # --- iterate over orientations & sides -----------------------------------
    for row, res in enumerate(res_sorted):
        orient = "original" if not res.orientation.is_flipped else "flipped"
            
        _axis(
            axs[row, 0],
            x=PHOSPHATE_SITES,
            tau=res.tau_left,
            dG=res.dG_left,
            side_label=f"(L, {orient})",
            color=colour[row][0],
        )
        axs[row, 0].set_title(f"τ & ΔG - left unwrap ({orient})")

        _axis(
            axs[row, 1],
            x=PHOSPHATE_SITES,
            tau=res.tau_right,
            dG=res.dG_right,
            side_label=f"(R, {orient})",
            color=colour[row][1],
        )
        axs[row, 1].set_title(f"τ & ΔG - right unwrap ({orient})")
        axs[row, 1].set_xlabel("Phosphate site (bp)")

    fig.suptitle(f"Dwell-time vs. step ΔG per orientation and side for sequence {seq_id}" , y=1.03)
    plt.tight_layout()
    plt.show()



######### Eads: CALCULATED ADSORPTION ENERGIES WITH FREE ENERGY(F), ENTHALPY AND NO ENERGY ######################


def plot_eads_elegant(results: List[EadsResult], *, seq_id:str, figsize: Tuple[int, int] = (20, 10)) -> None:
    """Elegant 2x3 comparison of basic/full/enthalpy energies.

    • Row_0 = original orientation (blue left, red right)
    • Row_1 = flipped   orientation (red  left, blue right)
    """

    def _single_axis(ax: plt.Axes, x: List[int], y_left: pd.Series, y_right: pd.Series, *, title: str, ylabel: str, label_left: str, label_right: str, colour_left: str, colour_right: str) -> None:
        ax.plot(x, y_left, "-o", label=label_left, linewidth=2, alpha=0.85, color=colour_left)
        ax.plot(x, y_right, "--x", label=label_right, linewidth=2, alpha=0.85, color=colour_right)
        ax.set_title(title)
        ax.set_xlabel("Phosphate site (bp)")
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle="--", alpha=0.6)
        sns.despine(ax=ax)




    if len(results) != 2:
        raise ValueError("Provide exactly two EadsResult objects (original & flipped)")

    # deterministic order: ORIGINAL first row, FLIPPED second
    res_sorted = sorted(results, key=lambda r: r.orientation.value)

    fig, axes = plt.subplots(2, 3, figsize=figsize, sharey=True)

    plot_map = [
        ("Eads_Base_left", "Eads_Base_right", r"$E_{ads}\,(k_B T)$", "Basic"),
        ("Eads_dG_left", "Eads_dG_right", r"$E_{ads}^{full}\,(k_B T)$", "Free‑energy"),
        ("Eads_dH_left", "Eads_dH_right", r"$E_{ads}^{enth}\,(k_B T)$", "Enthalpy"),
    ]

    colours = (
        ("tab:blue", "tab:red"),   # original row
        ("tab:red", "tab:blue"),   # flipped  row (swapped)
    )

    for row, res in enumerate(res_sorted):
        left_col, right_col = colours[row]
        orient_lbl = "Orig" if res.orientation is Orientation.ORIGINAL else "Flip"

        for col, (attr_L, attr_R, ylabel, ttl) in enumerate(plot_map):
            _single_axis(
                axes[row, col],
                PHOSPHATE_SITES,
                getattr(res, attr_L),
                getattr(res, attr_R),
                title=f"{ttl} ({orient_lbl})",
                ylabel=ylabel if col == 0 else "",
                label_left="left-unwrap",
                label_right="right-unwrap",
                colour_left=left_col,
                colour_right=right_col,
            )

    axes[0, 0].legend(frameon=True)
    fig.suptitle(f"Unzipping energy landscape for Sequence {seq_id}", y=1.02, fontsize=16)
    fig.tight_layout()
    plt.show()
