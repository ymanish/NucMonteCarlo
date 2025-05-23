import pandas as pd
from py_analysis.config.io_path import DATA_DIR
from py_analysis.adsorbpipe.utils.data_io import load_dwell_data
from py_analysis.config.eads_var import PHOSPHATE_SITES
from py_analysis.adsorbpipe.utils.helpers_plotting import plot_rawdwell_times, plot_cumulative_dwell_times


def calculate_cumulative_dwell_times(df: pd.DataFrame, 
                                    position_offset: int = 73,
                                    direction: str = 'left') -> pd.DataFrame:
    
    """Calculate cumulative dwell times with position adjustment."""
    df = df.copy()
    df["bp_pos"] = df["basepair"] + position_offset
    df.loc[df['dwelltime'] < 0, 'dwelltime'] = 0
    df = df.sort_values('bp_pos')
    
    if direction == 'right':
        df['cumulative_dwell_time'] = df['dwelltime'][::-1].cumsum()[::-1]
    else:
        df['cumulative_dwell_time'] = df['dwelltime'].cumsum()
    
    return df


def get_phosphate_sites(df_left: pd.DataFrame, 
                           df_right: pd.DataFrame) -> pd.DataFrame:
    
    """Process and merge phosphate site data."""
    # Process left zip data
    left_phos_df = df_left[df_left['bp_pos'].isin(PHOSPHATE_SITES)]
    
    ## Subtract 1 from the bp_pos because when zipping from right the bp will be of right to the bound phosphate. And since in the phosphate_sites 
    ## we are using the left base pair to represent the phosphate site, I need to make sure I select the dwell time for the base pair that is in the right 
    ## of the phosphate site, when I zip from the right.
    right_phos_df = df_right.copy()
    right_phos_df["aligned_bp_position"] = right_phos_df["bp_pos"] - 1
    right_phos_df = right_phos_df[right_phos_df['aligned_bp_position'].isin(PHOSPHATE_SITES)]
    
    merged_df = pd.merge(
        left_phos_df,
        right_phos_df,
        left_on='bp_pos',
        right_on='aligned_bp_position',
        suffixes=('_left', '_right')
    )
    
    final_df = merged_df[[
        "bp_pos_left", 
        "cumulative_dwell_time_left", 
        "cumulative_dwell_time_right"
    ]].rename(columns={
        "bp_pos_left": "phosphate_site",
        "cumulative_dwell_time_left": "left_cumdwell_time",
        "cumulative_dwell_time_right": "right_cumdwell_time"
    })
    
    return final_df

def getcumulative_dwelltimes(leftzip_filepath:str, right_filepath:str, plot:bool=False) -> pd.DataFrame:

    left_df, right_df = load_dwell_data(left_file=leftzip_filepath, 
                                        right_file=right_filepath)
    
    left_cum = calculate_cumulative_dwell_times(left_df, direction='left')
    right_cum = calculate_cumulative_dwell_times(right_df, direction='right')
        
    final_df = get_phosphate_sites(left_cum, right_cum)
    
    if plot:
        plot_rawdwell_times(left_df, right_df)
        plot_cumulative_dwell_times(left_cum, right_cum)

    return final_df

if __name__ == "__main__":
    
    left_path = DATA_DIR / f"Hall_forward_zip_dwelltime_clean.csv"
    right_path = DATA_DIR / f"Hall_reverse_zip_dwelltime_clean.csv"

    final_results = getcumulative_dwelltimes(leftzip_filepath=left_path,  
                                             right_filepath=right_path, 
                                             plot=True)
    print(final_results.head())