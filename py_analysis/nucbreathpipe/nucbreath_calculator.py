import sys
from py_analysis.modules.NucFreeEnergy import NucleosomeBreath
from py_analysis.config.custom_types import NuclBreathingResult
from typing import List, Tuple, Optional



def BreathEnergies_per_sequence(key, seq, 
                                nucmethod:str, 
                                bind_sates:List[Tuple[int, int]], 
                                factor:float,
                                hard:bool=False, 
                                style:str="b_index",
                                flip:bool=False, 
                                freedna_method: Optional[str] = None)->  List[NuclBreathingResult]:
    nucleosomebreath = NucleosomeBreath(nuc_method=nucmethod, free_dna_method=freedna_method)
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
