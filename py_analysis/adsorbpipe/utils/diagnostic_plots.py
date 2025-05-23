from py_analysis.config.custom_types import DwellTimeResult, FitResult
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, Sequence

from py_analysis.adsorbpipe.utils.helpers_plotting import plot_detected_peaks, plot_gaussian_fit, plot_shifted_peaks

def _create_peak_detection_fig(leftpeaks:Sequence[int],
                                rightpeaks:Sequence[int], 
                                org_left_dwelltime:np.ndarray | pd.Series,
                                org_right_dwelltime:np.ndarray | pd.Series,
                                bp_range:np.ndarray) -> None:
    """Figure 1: Raw data with detected peaks (left/right)"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))

    # Left zip
    plot_detected_peaks(
        dwelltime=org_left_dwelltime,
        bp_pos=bp_range,
        peaks=leftpeaks,
        title="Left Zip: Peak Detection",
        ax=ax1
    )
    
    # Right zip
    plot_detected_peaks(
        dwelltime=org_right_dwelltime,
        bp_pos=bp_range,
        peaks=rightpeaks,
        title="Right Zip: Peak Detection",
        ax=ax2
    )
    
    plt.tight_layout()
    plt.show()

def _create_gaussian_fit_fig(left_fit_result:FitResult, 
                             right_fit_result:FitResult, 
                             org_left_dwelltime:np.ndarray | pd.Series, 
                             org_right_dwelltime:np.ndarray | pd.Series, 
                             bp_range:np.ndarray):
    """Figure 2: Gaussian decomposition fits (left/right)"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))
    # Left zip
    plot_gaussian_fit(
        dwelltime=org_left_dwelltime,
        bp=bp_range,
        fit_result=left_fit_result,
        title="Left Zip: Gaussian Fit",
        ax=ax1
    )
    
    # Right zip
    plot_gaussian_fit(
        dwelltime=org_right_dwelltime,
        bp=bp_range,
        fit_result=right_fit_result,
        title="Right Zip: Gaussian Fit",
        ax=ax2
    )
    
    plt.tight_layout()
    plt.show()

def _create_peak_assignment_fig(leftpeaks:Sequence[int],
                                rightpeaks:Sequence[int], 
                                left_shited_peaks:dict[int, float],
                                right_shited_peaks:dict[int, float],
                                org_left_dwelltime:np.ndarray | pd.Series, ### Need for the original dwell time trace
                                org_right_dwelltime:np.ndarray | pd.Series,  ### Need for the original dwell time trace
                                bp_range:np.ndarray, 
                                left_prioroffset:float=0.0,
                                right_prioroffset:float=0.0) -> None:
    
    """Figure 3: Peak assignments (left/right)"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))

    # Left zip
    plot_shifted_peaks(
        dwelltime=org_left_dwelltime,
        bp=bp_range,
        detected_peaks=leftpeaks,
        shifted_sites=list(left_shited_peaks.keys()),
        values=list(left_shited_peaks.values()),
        offset=left_prioroffset,
        title="Left Zip: Peak Assignment",
        ax=ax1
    )
    
    # Right zip
    plot_shifted_peaks(
        dwelltime=org_right_dwelltime,
        bp=bp_range,
        detected_peaks=rightpeaks,
        shifted_sites=list(right_shited_peaks.keys()),
        values=list(right_shited_peaks.values()),
        offset=right_prioroffset,
        title="Right Zip: Peak Assignment",
        ax=ax2
    )
    
    plt.tight_layout()
    plt.show()


