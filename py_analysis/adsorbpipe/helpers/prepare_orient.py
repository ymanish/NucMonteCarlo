from typing import Tuple
import pandas as pd
from py_analysis.config.custom_types import Orientation


def _prepare_orientation(tau_left:pd.Series, tau_right:pd.Series, 
                         left_dG_df:pd.DataFrame, right_dG_df:pd.DataFrame, 
                         orientation: Orientation) -> Tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    """Return tau and deltaG/deltaH arrays mapped into *physical* left -> right order.
    tau_left: pd.Series  # from experimental unwrapping -73 -> +73
    tau_right: pd.Series  # from experimental unwrapping +73 -> -73 
    left_dG_df: pd.DataFrame  # Left is with respect to the histone octamer orientation example: 0 -> 27 
    right_dG_df: pd.DataFrame  # Right is with respect to the histone octamer orientation example: 0 -> 27
    orientation: Orientation  # physical orientation of the DNA on the histone octamer"""



      # ---------- LEFT SIDE ENERGY (RELATIVE TO DYAD)  ------------------
    dG_L  = left_dG_df['dF'].diff().shift(-1)         
    dH_L = left_dG_df['F_enthalpy'].diff().shift(-1)         

    dG_L.iloc[-1] = 0
    dH_L.iloc[-1] = 0  

    dG_L, dH_L  = dG_L.reset_index(drop=True), dH_L.reset_index(drop=True)                



    # ---------- RIGHT SIDE ENERGY (RELATIVE TO DYAD) ------------------
    tempF_R = right_dG_df['dF'][::-1].reset_index(drop=True)
    tempF_R_enthalpy = right_dG_df['F_enthalpy'][::-1].reset_index(drop=True)

    dG_R  = tempF_R.diff().shift(-1)                   
    dG_R.iloc[-1] = 0

    dH_R  = tempF_R_enthalpy.diff().shift(-1)                   
    dH_R.iloc[-1] = 0

    dG_R, dH_R  = dG_R[::-1].reset_index(drop=True), dH_R[::-1].reset_index(drop=True)

    if orientation.is_flipped:
          tau_L = tau_right[::-1].reset_index(drop=True)
          tau_R = tau_left[::-1].reset_index(drop=True)

    else:
          tau_L = tau_left.reset_index(drop=True)
          tau_R = tau_right.reset_index(drop=True)

    return tau_L, tau_R, dG_L, dG_R, dH_L, dH_R