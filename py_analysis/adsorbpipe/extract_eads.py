
from py_analysis.adsorbpipe.core.compute_Eads import compute_Eads_values, symmetrize_Eads
from py_analysis.config.io_path import DATA_DIR, RESULTS_DIR
from py_analysis.adsorbpipe.utils.data_io import load_dwell_data, load_dna_histone_energies
from py_analysis.adsorbpipe.core.gaussian_dwelltime import get_gaussian_dwelltimes
from py_analysis.config.eads_var import Parameters
from py_analysis.adsorbpipe.helpers.process_data import get_left_right_breathing, interpolate_dwelltime_to_bp
import numpy as np
from py_analysis.config.custom_types import Orientation




if __name__ == "__main__":


    SEQUENCE = "601" ### This is the only sequence available for unzipping dwell times
    params = Parameters(NUCMETHOD="crystal",
                        FREEDNA_METHOD="md")

    # left_path = DATA_DIR / f"Hall_forward_zip_dwelltime_clean.csv"
    # right_path = DATA_DIR / f"Hall_reverse_zip_dwelltime_clean.csv"

    left_path = DATA_DIR / f"left.csv"
    right_path = DATA_DIR / f"right.csv"



    left_df, right_df = load_dwell_data(left_path, right_path)
    left_df, right_df = interpolate_dwelltime_to_bp(left_df, right_df)


    ######## Get the dwell time ##########################

    results = get_gaussian_dwelltimes(left_df=left_df,
                                    right_df=right_df,
                                    peak_params = {
                                        "height": 1e-3,
                                        "prominence": 3e-4,
                                        "distance": 1},
                                    left_prioroffset = 0.0,
                                    right_prioroffset = -0.5) 

    dwell_times_df = results.to_dataframe()


    dna_histone_dir = RESULTS_DIR / f"nbfiles/dna_histone"
    df_full, df_full_flip = load_dna_histone_energies(dna_histone_dir=dna_histone_dir,
                                                        nucmethod=params.NUCMETHOD,
                                                        freedna_method=params.FREEDNA_METHOD,
                                                        flip=params.FLIP,
                                                        krescfactor=params.KRESCFACTOR)




    ### SELECT ONLY THE 601 SEQUENCE
    ### Because the unzipping dwell times are available only for the 601 sequence

    df_full = df_full[df_full["id"]==SEQUENCE]
    df_full_flip = df_full_flip[df_full_flip["id"]==SEQUENCE]

    #############################################
    #####ORIGINAL SEQUENCE ######################
    #############################################

    df_org_L, df_org_R = get_left_right_breathing(df_full)

    #############################################
    #####FLIPPED SEQUENCE #######################
    #############################################

    df_flip_L, df_flip_R = get_left_right_breathing(df_full_flip)    




    EadsResult_org = compute_Eads_values(Tau=dwell_times_df, 
                                forward_dG_df=df_org_L, 
                                reverse_dG_df=df_org_R, 
                                orientation=Orientation.ORIGINAL,
                                tau_0=params.TAU0)

    EadsResult_flipped = compute_Eads_values(Tau=dwell_times_df,
                                        forward_dG_df=df_flip_L, 
                                        reverse_dG_df=df_flip_R, 
                                        orientation=Orientation.FLIPPED,
                                        tau_0=params.TAU0)


    Eads_values = symmetrize_Eads(EadsResult_org, EadsResult_flipped, padding_edges=True)

    np.save(RESULTS_DIR / f"eads/param/nuc{params.NUCMETHOD}_freedna_{params.FREEDNA_METHOD}_K{params.KRESCFACTOR}_Eads_values.npy", Eads_values)