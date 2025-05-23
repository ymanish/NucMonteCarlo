from typing import Tuple
import pandas as pd
from pathlib import Path
import warnings

def load_dwell_data(left_file: str, right_file: str, position_offset:int=73) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess dwell time data from CSV files."""
    left_df = pd.read_csv(left_file, header=None, names=["basepair", "dwelltime"])
    right_df = pd.read_csv(right_file, header=None, names=["basepair", "dwelltime"])
    left_df["bp_pos"] = left_df["basepair"]+position_offset
    right_df["bp_pos"] = right_df["basepair"]+position_offset
    left_df.loc[left_df['dwelltime'] < 0, 'dwelltime'] = 0
    right_df.loc[right_df['dwelltime'] < 0, 'dwelltime'] = 0
    return left_df, right_df

def load_dna_histone_energies(dna_histone_dir:Path, 
                              nucmethod:str,
                              freedna_method: str, 
                              krescfactor:float, 
                              flip:bool = True,
                                hard_constraints: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # RESULTS_DIR / f"nbfiles/dna_histone

    if hard_constraints:
        ############## HARD CONSTRAITNS ######################  
        df_full = pd.read_csv(f"{dna_histone_dir}/nuc_{nucmethod}_freedna_{freedna_method}_breathstatefe_K{krescfactor}_hc.csv")
    else:
        ############## SOFT CONSTRAINTS ######################
        df_full = pd.read_csv(f"{dna_histone_dir}/nuc_{nucmethod}_freedna_{freedna_method}_breathstatefe_K{krescfactor}_sc.csv")

    if flip:
        #############################################
        #####FLIPPED SEQUENCE #######################
        #############################################
        df_full_flip = pd.read_csv(f"{dna_histone_dir}/nuc_{nucmethod}_freedna_{freedna_method}_breathstatefe_K{krescfactor}_sc_flip.csv")
    else:
        print("No flipped sequence data requested. Returning empty DataFrame.")
        warnings.warn(
            "No flipped sequence data requested; flipped data is important "
            "for adsorption energy calculations.",
            UserWarning
        )
        df_full_flip = pd.DataFrame()

    return df_full, df_full_flip



if __name__ == "__main__":
    ## Test the functions
    
    from py_analysis.config.io_path import DATA_DIR, RESULTS_DIR

            
    left_path = DATA_DIR / f"Hall_forward_zip_dwelltime_clean.csv"
    right_path = DATA_DIR / f"Hall_reverse_zip_dwelltime_clean.csv"
    left_df, right_df = load_dwell_data(left_path, right_path)

    dna_histone_dir = RESULTS_DIR / f"nbfiles/dna_histone"
    df_full, df_full_flip = load_dna_histone_energies(dna_histone_dir=dna_histone_dir,
                                                      nucmethod="hybrid",
                                                      freedna_method=None,
                                                      flip=False,
                                                      krescfactor=1.0)