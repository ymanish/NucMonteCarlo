from typing import Tuple
import warnings
import pandas as pd
import numpy as np
from py_analysis.config.custom_types import Orientation, EadsResult
from py_analysis.adsorbpipe.helpers.prepare_orient import _prepare_orientation
from py_analysis.config.eads_var import PHOSPHATE_SITES

def single_eads(tau:pd.Series,
                 dG:pd.Series, 
                 dH:pd.Series, 
                 kBT:float, 
                 tau_0:float) -> Tuple[float, float]:
    
    Eads_base = kBT*np.log(tau / tau_0) 
    Eads_dG = Eads_base - dG 
    Eads_dH = Eads_base- dH 
    return Eads_base, Eads_dG, Eads_dH


def compute_Eads_values(Tau:pd.DataFrame, 
                    forward_dG_df:pd.DataFrame,
                      reverse_dG_df:pd.DataFrame,
                       orientation:Orientation, 
                         tau_0:float=1e-6, kBT:float=1.0)-> EadsResult:
    
    # ---------- LEFT  (forward zipping) From Experiment on the original orientation of 601 ------------------
    tau_L = Tau["left_dwell_time"].reset_index(drop=True)               
   
    # ---------- RIGHT  (reverse zipping) From Experiment on the original orientation of 601 ------------------
    tau_R = Tau["right_dwell_time"].reset_index(drop=True)               

    

    tau_L, tau_R, dG_L, dG_R, dH_L, dH_R = _prepare_orientation(tau_left=tau_L, tau_right=tau_R,
                                                                                     left_dG_df=forward_dG_df, right_dG_df=reverse_dG_df,
                                                                                     orientation=orientation)

    Eads_L, Eads_L_dG, Eads_L_dH = single_eads(tau=tau_L, dG=dG_L, dH=dH_L, kBT=kBT, tau_0=tau_0)
    Eads_R, Eads_R_dG, Eads_R_dH = single_eads(tau=tau_R, dG=dG_R, dH=dH_R, kBT=kBT, tau_0=tau_0)
    

    # sanitise infinities caused by zero dwell times
    for s in (Eads_L, Eads_R, Eads_L_dG, Eads_R_dG, Eads_L_dH, Eads_R_dH):
        s.replace([np.inf, -np.inf], 0, inplace=True)


    return EadsResult(
        orientation,
        Eads_Base_left=Eads_L,
        Eads_Base_right=Eads_R,
        Eads_dG_left=Eads_L_dG,
        Eads_dG_right=Eads_R_dG,
        Eads_dH_left=Eads_L_dH,
        Eads_dH_right=Eads_R_dH,
        tau_left=tau_L,
        tau_right=tau_R,
        dG_left=dG_L,
        dG_right=dG_R,
        dH_left=dH_L,
        dH_right=dH_R
    )



def symmetrize_Eads(EadsResult_org: EadsResult, EadsResult_flipped: EadsResult, padding_edges:bool=True) -> np.ndarray:
    """
    Symmetrize the Eads results using the original and flipped orientations Eads data.
    """

    BPS_PER_SITE = 2                    # two bp per histone–DNA contact
    HALF = len(PHOSPHATE_SITES) // 2   

    def _reduce_to_14_sites(Eads: pd.DataFrame) -> pd.Series:
        """
        Reduce the Eads data from 28 phosphate sites to 14 binding sites.
        """
        # ------------------------------------------------------------------
        # keep only the *relevant* half for each side
        # (left curves describe sites 0...13, right curves 14...27)
        # ------------------------------------------------------------------
        left_mask = np.r_[np.ones(HALF, bool), np.zeros(HALF, bool)]
        right_mask = ~left_mask


        Eads_df_all.loc[right_mask, ["org_left_Eads",  "flip_left_Eads_(org_right_Eads_flipped)"]]  = 0
        Eads_df_all.loc[left_mask,  ["org_right_Eads", "flip_right_Eads_(org_left_Eads_flipped)"]] = 0
        Eads_mean = Eads_df_all.mean(axis=1)

        # ------------------------------------------------------------------
        # aggregate to *binding-site* resolution (2 bp per site)
        # ------------------------------------------------------------------
        Eads_agg = (
            Eads_mean
            .groupby(Eads_mean.index // BPS_PER_SITE)
            .sum()
        )

        return Eads_agg


    Eads_df_all = pd.concat([
        EadsResult_org.Eads_dG_left,
        EadsResult_org.Eads_dG_right,
        EadsResult_flipped.Eads_dG_left,
        EadsResult_flipped.Eads_dG_right
    ], axis=1, ignore_index=True)
    Eads_df_all.columns = ["org_left_Eads", "org_right_Eads", "flip_left_Eads_(org_right_Eads_flipped)", "flip_right_Eads_(org_left_Eads_flipped)"]

    Eads_bnd_site = _reduce_to_14_sites(Eads_df_all)

    if padding_edges:
    # edge padding: copy the third site into sites 0–1
        # and the fourth-from-last into the last two
        Eads_bnd_site.iloc[:2]  = Eads_bnd_site.iloc[3]
        Eads_bnd_site.iloc[-2:] = Eads_bnd_site.iloc[-4]
    
    return np.array(Eads_bnd_site.values) 

###############DEPRECATED######################
###############################################
###############################################
###############################################


def compute_Eads_from_cumtau(Tau_cum_df:pd.DataFrame, 
                            forward_dG_df:pd.DataFrame,
                            reverse_dG_df:pd.DataFrame,
                            orientation:Orientation, 
                                tau_0:float=1e-6, kBT:float=1.0)-> EadsResult:
    """
    [Deprecated] Calculate Eads values from cumulative dwell times.
    
    .. deprecated:: 1.0
        Use 'compute_Eads_values' instead. This method will be removed in future versions.
    """
    warnings.warn(
        "'get_Eads_values' is deprecated and will be removed in a future version. "
        "Use 'compute_Eads_values' instead.",
        DeprecationWarning,
        stacklevel=2
    )


    # ---------- LEFT  (forward zipping) From Experiment on the original orientation of 601 ------------------
    tau_L = Tau_cum_df['left_cumdwell_time'].diff().fillna(Tau_cum_df['left_cumdwell_time'].iloc[0]) # first row has no upstream neighbour, so filled with its own value
    tau_L = tau_L.reset_index(drop=True)               
   


    # ---------- RIGHT  (reverse zipping) From Experiment on the original orientation of 601 ------------------
    cumT_R = Tau_cum_df['right_cumdwell_time'][::-1].reset_index(drop=True)
    tau_R = cumT_R.diff().fillna(cumT_R.iloc[0])      
    tau_R = tau_R[::-1].reset_index(drop=True)
    

    tau_L, tau_R, dG_L, dG_R, dH_L, dH_R = _prepare_orientation(tau_left=tau_L, tau_right=tau_R,
                                                                                     left_dG_df=forward_dG_df, right_dG_df=reverse_dG_df,
                                                                                     orientation=orientation)

    Eads_L, Eads_L_dG, Eads_L_dH = single_eads(tau=tau_L, dG=dG_L, dH=dH_L, kBT=kBT, tau_0=tau_0)
    Eads_R, Eads_R_dG, Eads_R_dH = single_eads(tau=tau_R, dG=dG_R, dH=dH_R, kBT=kBT, tau_0=tau_0)
    

    # sanitise infinities caused by zero dwell times
    for s in (Eads_L, Eads_R, Eads_L_dG, Eads_R_dG, Eads_L_dH, Eads_R_dH):
        s.replace([np.inf, -np.inf], 0, inplace=True)


    return EadsResult(
        orientation,
        Eads_Base_left=Eads_L,
        Eads_Base_right=Eads_R,
        Eads_dG_left=Eads_L_dG,
        Eads_dG_right=Eads_R_dG,
        Eads_dH_left=Eads_L_dH,
        Eads_dH_right=Eads_R_dH,
        tau_left=tau_L,
        tau_right=tau_R,
        dG_left=dG_L,
        dG_right=dG_R,
        dH_left=dH_L,
        dH_right=dH_R
    )

