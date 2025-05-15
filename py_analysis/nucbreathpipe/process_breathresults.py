import pandas as pd
from py_analysis.config.breath_var import Parameters
from py_analysis.config.io_path import RESULTS_DIR
from py_analysis.config.custom_types import NuclBreathingResult
from typing import List

def process_results(results_all: List[NuclBreathingResult], params: Parameters) -> pd.DataFrame:

    df_results = pd.DataFrame([result._asdict() for result in results_all])
    df_free_energy = df_results['F_vals'].apply(lambda x: x._asdict() if hasattr(x, '_asdict') else x).apply(pd.Series)
    df_free_energy = df_free_energy[['F', 'F_entropy', 'F_enthalpy', 'F_freedna']]

    df_full = pd.concat([df_results.drop(columns=['F_vals']), df_free_energy], axis=1)
    df_full["dF"] = df_full["F"] - df_full["F_freedna"]


    if params.STYLE in ["b_index", "ph_index"]:
        df_full["left_open"] = df_full["leftbind_indx"]
        df_full["right_open"] = df_full["rightbind_indx"].apply(lambda x: (params.LENGTH - 1) - x)
    else:
        df_full["left_open"] = df_full["leftbind_indx"]
        df_full["right_open"] = df_full["rightbind_indx"]

        
    return df_full

def save_results(df_full: pd.DataFrame, params: Parameters, dir:str="nbfiles") -> None:
    if params.FOR_DNA_HISTONE:
        subdir = "dna_histone"
    else:
        subdir = "nucbreathfe"

    if params.nucmethod == "hybrid":
        method_str = "hybrid"
    elif params.nucmethod == "crystal":
        method_str = "crystal"
    elif params.nucmethod == "md":
        method_str = "md"
    else:
        raise ValueError(f"Unknown nucleosome method: {params.nucmethod}")

    if params.HARD_CONS:
        filename = f"{method_str}_breathstatefe_K{params.KRESCFACTOR}_hc.csv"
    else:
        flip_str = "_flip" if params.FLIP_SEQUENCE else ""
        filename = f"{method_str}_breathstatefe_K{params.KRESCFACTOR}_sc{flip_str}.csv"


    filepath = RESULTS_DIR / dir / subdir / filename
    df_full.to_csv(filepath, index=False)