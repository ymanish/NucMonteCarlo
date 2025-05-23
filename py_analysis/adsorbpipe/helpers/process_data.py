from py_analysis.config.eads_var import PHOSPHATE_SITES
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

def get_left_right_breathing(df):

    """Separates the breathing purely from left and right side and adds the phosphate site information"""

    df_left = df[df['right_open']==0].reset_index(drop=True)
    df_right = df[df['left_open']==0].reset_index(drop=True)

    df_left['phosphate_site'] = np.where(
        df_left['leftbind_indx'] < len(PHOSPHATE_SITES),
        df_left['leftbind_indx'].astype(int).map(lambda i: PHOSPHATE_SITES[i]),
        np.nan
    )

    df_right['phosphate_site'] = np.where(
        df_right['rightbind_indx'] < len(PHOSPHATE_SITES),
        df_right['rightbind_indx'].astype(int).map(lambda i: PHOSPHATE_SITES[i]),
        np.nan
    )

    df_left = df_left[["phosphate_site", "id", "F", "F_entropy", "F_enthalpy", "F_freedna", "dF"]]
    df_right = df_right[["phosphate_site", "id", "F", "F_entropy", "F_enthalpy", "F_freedna", "dF"]]

    return df_left, df_right


def interpolate_dwelltime_to_bp(left_df, right_df):

    # Interpolate to base pair resolution
    bps = np.arange(np.ceil(min(left_df['bp_pos'].min(), right_df['bp_pos'].min())),
                    np.floor(max(left_df['bp_pos'].max(), right_df['bp_pos'].max())) + 1, 1)

    interp_left = interp1d(left_df['bp_pos'], left_df['dwelltime'], kind='linear', fill_value="extrapolate")
    interp_right = interp1d(right_df['bp_pos'], right_df['dwelltime'], kind='linear', fill_value="extrapolate")

    # Combine interpolated data
    combined_df = pd.DataFrame({
        'bp_pos': bps,
        'dwell_left': interp_left(bps),
        'dwell_right': interp_right(bps),
        'dwell_total': interp_left(bps) + interp_right(bps)
    })

    left_df = combined_df[['bp_pos', 'dwell_left']]
    right_df = combined_df[['bp_pos', 'dwell_right']]
    left_df.columns = ['bp_pos', 'dwelltime']
    right_df.columns = ['bp_pos', 'dwelltime']
    left_df["bp_pos"] = left_df["bp_pos"].astype(int)
    right_df["bp_pos"] = right_df["bp_pos"].astype(int)

    return left_df, right_df