import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from typing import Sequence, Tuple
from py_analysis.config.eads_var import PHOSPHATE_SITES

def detect_peaks(dwelltime: np.ndarray | pd.Series,
                    bp_pos: np.ndarray,
                    *,
                    height: float = 1e-2,
                    prominence: float = 3e-4,
                    distance: int = 1) -> np.ndarray:
    
    """Return positions (in bp) of detected peaks."""
    indices, _ = find_peaks(dwelltime,
                            height=height,
                            prominence=prominence,
                            distance=distance)
    return bp_pos[indices]

def assign_peaks_to_sites(peak_table: pd.DataFrame, priorpeak_offset:float=0.0) -> Tuple[dict[int, float], dict[int, float]]:
    ph_sites = np.array(PHOSPHATE_SITES)
    """Return {site : dwell_time} ensuring each site is used at most once."""
    # avail: set[int] = set(ph_sites)
    dwell: dict[int, float] = {s: 0.0 for s in ph_sites}
    shifted_peaks: dict[int, float] = {}## for plotting later to see which peak moved to which phosphate site

    for _, row in peak_table.iterrows():
        peak_pos = float(row.peak_pos) + priorpeak_offset
        distances = np.abs(ph_sites - peak_pos)
        nearest_idx = np.argmin(distances)
        nearest_site = ph_sites[nearest_idx]

        dwell[nearest_site] = row.dwell_time
        shifted_peaks[nearest_site] = row.A ### for plotting: "A" is the height of the peak

    return dwell, shifted_peaks


