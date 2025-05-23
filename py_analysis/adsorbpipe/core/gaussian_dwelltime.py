from __future__ import annotations

from py_analysis.adsorbpipe.utils.data_io import load_dwell_data
from py_analysis.adsorbpipe.helpers.peaks_fnc import detect_peaks, assign_peaks_to_sites
from py_analysis.adsorbpipe.helpers.fit_multigaussian import fit_multi_gaussian
from py_analysis.config.eads_var import PHOSPHATE_SITES
from py_analysis.config.io_path import DATA_DIR
from py_analysis.config.custom_types import DwellTimeResult
import numpy as np
import pandas as pd


def get_gaussian_dwelltimes(left_df: pd.DataFrame,
                            right_df: pd.DataFrame,
                            peak_params: dict = {
                                "height": 1e-3,
                                "prominence": 3e-4,
                                "distance": 1
                                                },
                            left_prioroffset:float = 1.0,
                            right_prioroffset:float = -1.0) -> DwellTimeResult:
    
    """Perform Gaussian fitting and peak detection on raw dwell time data."""
        
    bp_range = left_df["bp_pos"].values
    # validate that the two arrays are identical
    if not np.array_equal(bp_range, right_df["bp_pos"].values):
        raise ValueError("Left and right dataframes must have the same base pair positions.")

    org_left_dwelltime = left_df.dwelltime.values
    org_right_dwelltime = right_df.dwelltime.values

    
    left_peaks = detect_peaks(org_left_dwelltime, bp_range, **peak_params)
    right_peaks = detect_peaks( org_right_dwelltime, bp_range, **peak_params)
    print(f"Left peaks: {left_peaks}")
    print(f"Right peaks: {right_peaks}")

    left_fit = fit_multi_gaussian(org_left_dwelltime, bp_range, left_peaks)
    right_fit = fit_multi_gaussian(org_right_dwelltime, bp_range, right_peaks)
    
    # Map the dwell times calculated to nearest phosphate sites
    left_mapped, shift_left_A_ = assign_peaks_to_sites(left_fit.results, priorpeak_offset=left_prioroffset)
    right_mapped, shift_right_A_ = assign_peaks_to_sites(right_fit.results, priorpeak_offset=right_prioroffset)
    
    # Create final container
    sites = sorted(PHOSPHATE_SITES)

    return DwellTimeResult(
        phosphate_sites=np.array(sites),
        left_dwell_times=np.array([left_mapped.get(s, 0.0) for s in sites]),
        right_dwell_times=np.array([right_mapped.get(s, 0.0) for s in sites]),
    )




if __name__ == "__main__":

    from py_analysis.adsorbpipe.utils.diagnostic_plots import _create_gaussian_fit_fig, _create_peak_detection_fig, _create_peak_assignment_fig

        
    left_path = DATA_DIR / f"Hall_forward_zip_dwelltime_clean.csv"
    right_path = DATA_DIR / f"Hall_reverse_zip_dwelltime_clean.csv"

    left_df, right_df = load_dwell_data(left_path, right_path)

    bp_range = left_df["bp_pos"].values
    peak_params = {
        "height": 1e-3,
        "prominence": 3e-4,
        "distance": 1
    }
    left_prioroffset=+1.0
    right_prioroffset=-1.0
    
    org_left_dwelltime = left_df.dwelltime.values
    org_right_dwelltime = right_df.dwelltime.values

    
    left_peaks = detect_peaks(org_left_dwelltime, bp_range, **peak_params)
    right_peaks = detect_peaks( org_right_dwelltime, bp_range, **peak_params)
    print(f"Left peaks: {left_peaks}")
    print(f"Right peaks: {right_peaks}")

    left_fit = fit_multi_gaussian(org_left_dwelltime, bp_range, left_peaks)
    right_fit = fit_multi_gaussian(org_right_dwelltime, bp_range, right_peaks)

    print(f"Left fit results: {left_fit.results}")
    print(f"Right fit results: {right_fit.results}")

    
     # Map the dwell times calculated to nearest phosphate sites
    left_mapped, shift_left_A_ = assign_peaks_to_sites(left_fit.results, priorpeak_offset=left_prioroffset)
    right_mapped, shift_right_A_ = assign_peaks_to_sites(right_fit.results, priorpeak_offset=right_prioroffset)


    # Generate plots (optional)


    # _create_peak_detection_fig(leftpeaks=left_peaks,
    #                             rightpeaks=right_peaks,
    #                             org_left_dwelltime=left_df.dwelltime.values,
    #                             org_right_dwelltime=right_df.dwelltime.values, 
    #                              bp_range=bp_range)

    # _create_gaussian_fit_fig(left_fit_result=left_fit,
    #                         right_fit_result=right_fit, 
    #                         org_left_dwelltime=left_df.dwelltime.values,
    #                         org_right_dwelltime=right_df.dwelltime.values, 
    #                         bp_range=bp_range)
    
    _create_peak_assignment_fig(leftpeaks=left_peaks,
                                rightpeaks=right_peaks,
                                left_shited_peaks=shift_left_A_,
                                right_shited_peaks=shift_right_A_,
                                org_left_dwelltime=left_df.dwelltime.values,
                                org_right_dwelltime=right_df.dwelltime.values, 
                                bp_range=bp_range, 
                                left_prioroffset=left_prioroffset,
                                right_prioroffset=right_prioroffset)


        

    # results = get_gaussian_dwelltimes(left_df, right_df)

    # final_df = results.to_dataframe()
    # print(final_df)