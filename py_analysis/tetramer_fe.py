

from py_analysis.modules.NucFreeEnergy import NucleosomeBreath ## Needs to be import before numpy and numba because modules.NucFreeEnergy import env_settings 
#which sets the global limit of multithreading for each process.

import numpy as np
import numba as nb
import time
import concurrent.futures
from tqdm import tqdm
from typing import List, Tuple
import os
import sys
import pickle


from py_analysis.config.custom_types import ProcessedSequence, FreeEnergyResult
from py_analysis.config.seq_var import *
from py_analysis.config.gen_var import *

def process_sequence_padding_and_sliding(seq:str, 
                                         seq_id:str, 
                                         window_length: int, 
                                         final_length: int = 147, 
                                         pad_char: str = "A")-> List[ProcessedSequence]:
    
    if window_length > final_length:
        raise ValueError("Window length cannot be greater than the final desired length.")

    total_pad = final_length - window_length
    left_pad = total_pad // 2
    right_pad = total_pad - left_pad

    windows: List[ProcessedSequence] = []
    for i in range(len(seq) - window_length + 1):
        window = seq[i:i + window_length]
        padded_window = pad_char * left_pad + window + pad_char * right_pad
        # subid can be defined as sequence_id plus the window start index
        subid = f"{seq_id}_{i}"
        processed_seq = ProcessedSequence(
            id=seq_id,
            subid=subid,
            sequence=padded_window,
            start_site=i,
            end_site=i + window_length - 1  # original window boundaries
        )
        windows.append(processed_seq)
    
    return windows



def energy_per_long_sequence(key, records, nucmethod:str, hard:bool=False, k_factor:float=1.0)-> Tuple[str, List[FreeEnergyResult]]:
    nucleosomebreath = NucleosomeBreath(nuc_method=nucmethod)
    tetramer_loc = BIND_POINTS
    results:List[FreeEnergyResult] = []
    for rec in records:

        if hard:
            free_energy = nucleosomebreath.calculate_free_energy_hard(seq147=rec.sequence,
                                                                    left=tetramer_loc[0], right=tetramer_loc[1], id=key, subid=rec.subid)

        else:
            free_energy = nucleosomebreath.calculate_free_energy_soft(seq601=rec.sequence, 
                                                                  left=tetramer_loc[0], right=tetramer_loc[1], id=key, subid=rec.subid,  kresc_factor=k_factor)
        results.append(free_energy)

    return key, results

# def cgnaplus_fn(seq):
#     gs,stiff = cgnaplus_bps_params(seq,group_split=True)
#     return gs,stiff

if __name__ == "__main__":
  
    start = time.perf_counter()

    KRESCFACTOR = float(sys.argv[1])

    processed_sequences = {}
    for name, seq in sequence_dict.items():
        processed_sequences[name] = process_sequence_padding_and_sliding(seq=seq, 
                                                                     seq_id=name,
                                                                     window_length=TETRAMER_LENGTH, 
                                                                     final_length=NUC_LENGTH, 
                                                                     pad_char=PAD_CHAR)
        




    results_all = dict()
    with concurrent.futures.ProcessPoolExecutor(max_workers=11) as executor:
        futures = [executor.submit(energy_per_long_sequence, key=rec, 
                                                            records=processed_sequences[rec], 
                                                            nucmethod=NUC_PARAM_TYPE, 
                                                            hard=HARD_CONS, 
                                                            k_factor=KRESCFACTOR) for rec in processed_sequences.keys()]
        total = len(futures)

        for future in tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences"):
            seq_id, energy_array = future.result()
            print(f"Appending results for sequence {seq_id}:")

            results_all[seq_id] = energy_array

    print (f"Total sequences processed: {len(results_all)}")
    # print (f"Results: {results}")


    if HARD_CONS:
        file_cfg = ddGDataSaveParams(
                        nuc_type=NUC_PARAM_TYPE,
                        hangdna_type=HANG_PARAM_TYPE,
                        tetramer_length=TETRAMER_LENGTH,
                        pad_char=PAD_CHAR,
                        constraint="hc", 
                        kresc_factor=KRESCFACTOR
                    )
       
    else:
        file_cfg = ddGDataSaveParams(
                nuc_type=NUC_PARAM_TYPE,
                hangdna_type=HANG_PARAM_TYPE,
                tetramer_length=TETRAMER_LENGTH,
                pad_char=PAD_CHAR,
                constraint="sc", 
                kresc_factor=KRESCFACTOR

            )

    pkl_filepath = RESULTS_DIR / "pklfiles" / f"{file_cfg.ddG_filename(long_name=False)}.pkl"
    with open(pkl_filepath, 'wb') as f:
        pickle.dump(results_all, f)

    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")