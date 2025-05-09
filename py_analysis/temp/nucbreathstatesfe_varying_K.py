import os 
import sys
from py_analysis.modules.NucFreeEnergy import NucleosomeBreath
from py_analysis.config.custom_types import FreeEnergyResult, NuclBreathingResult
from py_analysis.config.gen_var import RESULTS_DIR, HARD_CONS
from typing import List
import concurrent.futures
from tqdm import tqdm
import pandas as pd

def energy_per_sequence(key, seq, nucmethod:str, bind_sates:List[tuple], factor:float, hard:bool=False, style:str="b_index", flip:bool=False)->  List[NuclBreathingResult]:
    nucleosomebreath = NucleosomeBreath(nuc_method=nucmethod)
    results:List[NuclBreathingResult] = []
    print(seq)
    if flip:
        translation_table = str.maketrans("ATCG", "TAGC")
        seq = seq[::-1]
        seq = seq.translate(translation_table)
        
    print(f"Sequence: {seq}")

    for bind_loc in bind_sates:

        if hard:
            free_energy = nucleosomebreath.calculate_free_energy_hard(seq147=seq,
                                                                    left=bind_loc[0], right=bind_loc[1], id=key)

        else:
            free_energy = nucleosomebreath.calculate_free_energy_soft(seq601=seq, 
                                                                  left=bind_loc[0], right=bind_loc[1], id=key, kresc_factor=factor, style=style)
            
        res_nucbreath = NuclBreathingResult(
            id=key,
            subid=None,
            sequence=seq,
            leftbind_indx=bind_loc[0],
            rightbind_indx=bind_loc[1],
            F_vals=free_energy
        )
        results.append(res_nucbreath)

    return results

def total_states_index(length:int=14):
    """
    Returns the total number of states.
    """
    states = []
    for left in range(length):  # Left binding site: 0 to 13
        for right in range(left, length):  # Right must be ≥ left
            states.append((left, right))

    print(f"Total states index wise: {len(states)}")
    ## Use this for style="bi"
    return states


def total_open_states(length:int=28):
    """
    Returns the total number of open states.
    """
   
    open_states = []
    for left in range(29): 
        for right in range(29):  
            if left+right <= 28: # Left and right binding sites must not exceed 28
                open_states.append((left, right))


    print(f"Total open states: {len(open_states)}")
    ## Use this for style="os"
    return open_states


if __name__ == "__main__":


    #############################################################
    #####################PARAMETERS #############################
    #############################################################
    #############################################################
    

    seq_dict = {"601":"CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT", 
                "601RTA":"CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT", 
                "601MF": "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACACCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT",
                "601L": "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAG",
                "5S": "CTTCCAGGGATTTATAAGCCGATGACGTCATAACATCCCTGACCCTTTAAATAGCTTAACTTTCATCAAGCAAGAGCCTACGACCATACCATGCTGAATATACCGGTTCTCGTCCGATCACCGAAGTCAAGCAGCATAGGGCTCGGT"  }

    nucmethod = "hybrid"
    STYLE = "ph_index" # "b_index" or "open_sites" or "ph_index"
    FOR_DNA_HISTONE = True
    FLIP_SEQUENCE = True

    if STYLE == "open_sites":
        states = total_open_states()
    elif STYLE == "b_index":
        LENGTH = 14
        states = total_states_index(length=14)
    else:
        LENGTH = 28
        states = total_states_index(length=28)
    
    KRESCFACTOR = float(sys.argv[1])

    print(f"PARAMETERS...............")
    print(f"nucmethod: {nucmethod}")
    print(f"FOR DNA-HISTONE: {FOR_DNA_HISTONE}")
    print(f"Sequence Flipped: {FLIP_SEQUENCE}")
    print(f"Style: {STYLE}")
    print(f"Length: {LENGTH}")
    print(f"Factor: {KRESCFACTOR}")
    print(f"Hard constraints: {HARD_CONS}")
    print(f"Total states: {len(states)}")
    print(f"Total sequences: {len(seq_dict)}")
    #############################################################
    #####################MAIN BLOCK #############################
    #############################################################
    #############################################################

    results_all:List[NuclBreathingResult] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:


        futures = [executor.submit(energy_per_sequence, key=s, seq=seq_dict[s],
                                    nucmethod='hybrid', 
                                    bind_sates=states, 
                                    factor=KRESCFACTOR, 
                                    hard=HARD_CONS, style=STYLE, flip=FLIP_SEQUENCE) for s in seq_dict.keys()]        
        total = len(futures)

        for future in tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences"):
            energy_array = future.result()
            results_all.extend(energy_array)



    df_results = pd.DataFrame([result._asdict() for result in results_all])
    df_free_energy = df_results['F_vals'].apply(lambda x: x._asdict() if hasattr(x, '_asdict') else x).apply(pd.Series)
    df_free_energy = df_free_energy[['F', 'F_entropy', 'F_enthalpy', 'F_freedna']]
    df_full = pd.concat([df_results.drop(columns=['F_vals']), df_free_energy], axis=1)


    df_full["dF"] = df_full["F"] - df_full["F_freedna"]



    if STYLE == "b_index" or STYLE == "ph_index":
        df_full["left_open"] = df_full["leftbind_indx"]
        df_full["right_open"] = df_full["rightbind_indx"].apply(lambda x: (LENGTH-1) - x)

    else:
        df_full["left_open"] = df_full["leftbind_indx"]
        df_full["right_open"] = df_full["rightbind_indx"]



    if FOR_DNA_HISTONE:
        if HARD_CONS:
            df_full.to_csv(RESULTS_DIR / f"nbfiles/dna_histone/breathstatefe_K{KRESCFACTOR}_hc.csv", index=False)
        else:
            df_full.to_csv(RESULTS_DIR / f"nbfiles/dna_histone/breathstatefe_K{KRESCFACTOR}_sc_flip.csv", index=False)


    else:
        if HARD_CONS:
            df_full.to_csv(RESULTS_DIR / f"nbfiles/nucbreathfe/breathstatefe_K{KRESCFACTOR}_hc.csv", index=False)
        else:
            df_full.to_csv(RESULTS_DIR / f"nbfiles/nucbreathfe/breathstatefe_K{KRESCFACTOR}_sc.csv", index=False)

  

  