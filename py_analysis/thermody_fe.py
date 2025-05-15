from py_analysis.modules.NucFreeEnergy import NucleosomeBreath

import numba as nb
import numpy as np
from py_analysis.config.gen_var import TI_PARAMS
from py_analysis.config.io_path import DATA_DIR, RESULTS_DIR
import time 
import concurrent.futures
from tqdm import tqdm
import os
import pandas as pd
from typing import List, Optional
import sys
from py_analysis.config.custom_types import FreeEnergyResult

def process_batch(batch:List[dict], method_nuc:str)-> List[FreeEnergyResult]:
    
    nuc_breath = NucleosomeBreath(nuc_method=method_nuc)
    results:List[FreeEnergyResult]=[]

    # print(batch)
    for key,value in batch.items():
        # print("Processing sequence:", key, value)
        id = key
        subseq= value
        # F, F_entropy, E, F_free, F_diff = nuc_breath.calculate_free_energy(seq601=subseq, left_open=0, right_open=0)
        res = nuc_breath.calculate_free_energy_soft(seq601=subseq, left=0, right=13, id=id)
        results.append(res)
    return results


def batch_dict(data: dict, batch_size: int=10) -> List[dict]:

    items = list(data.items())
    return [dict(items[i:i+batch_size]) for i in range(0, len(items), batch_size)]


## process a reference sequence seqA 
def process_ref_seq(seqA:str, method_nuc, left:int=0, right:int=13, id:Optional[str]=None, subid:Optional[str]=None)-> List[FreeEnergyResult]:
    nuc_breath = NucleosomeBreath(nuc_method=method_nuc)
    free_energy = nuc_breath.calculate_free_energy_soft(seq601=seqA, left=left, right=right, id=id, subid=subid)
    return [free_energy]

if __name__ == "__main__":
  
    start = time.perf_counter()

    WHOLE_SEQ = TI_PARAMS['WHOLE_SEQ']
    SEQ_LENGTH = TI_PARAMS['SEQ_LENGTH']
    SLIDE_INTERVAL = TI_PARAMS['SLIDE_INTERVAL']
    MU_RANGE = TI_PARAMS['MU_RANGE']
    MU_VALUES = TI_PARAMS['MU_VALUES']
    results_all_batches:List[FreeEnergyResult]=[]

    if WHOLE_SEQ:
        parameter_file_path = DATA_DIR / f"S288C_YAL002W-flanking1k_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt"
    else:
        parameter_file_path = DATA_DIR / f"S288C_YAL002W-TSS_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt"



    # Read the parameter file (columns: ID, μ, SeqB_pos, SeqA, SeqB)
    df_param = pd.read_csv(parameter_file_path, sep='\t', header=None, names=['ID', 'mu', 'seqB_pos', 'seqA', 'seqB'])

    ###################### Process the reference sequence
    ref_seq = df_param['seqA'].iloc[0]    
    print (f"Reference sequence: {ref_seq}")

    # ref_results = process_ref_seq(ref_seq, method_nuc='hybrid', left=0, right=13, id='0', subid='ref_seq')

    # results_all_batches.extend(ref_results)
    # print (f"Reference sequence processed..")
    ##################### Process the batches of sequences
    nuc_breath_df = df_param[['seqB', 'seqB_pos']].drop_duplicates(subset=['seqB'], keep='first')
    seq_dict = {str(pos): seq for pos, seq in zip(nuc_breath_df['seqB_pos'], nuc_breath_df['seqB'])}
    batches = batch_dict(seq_dict, batch_size=100)


    with concurrent.futures.ProcessPoolExecutor(max_workers=11) as executor:

        # Process reference sequence as first batch
        # ref_batch = {"ref_seq": ref_seq}
        futures = [executor.submit(process_ref_seq, ref_seq, 'hybrid', 0, 13, '0', 'ref_seq')] 
        
        futures += [executor.submit(process_batch, batch, 'hybrid') for batch in batches]
        total = len(futures)

        for future in tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences"):
            energy_array = future.result()
            results_all_batches.extend(energy_array)
           
    print (f"Total sequences processed: {len(results_all_batches)}")
    # for i, result in enumerate(results_all_batches):
    #     print(f"Result {i}: {result}")
    results_dict_list = [r._asdict() for r in results_all_batches]

    # Create a DataFrame from the list of dictionaries
    results_df = pd.DataFrame(results_dict_list)


    if WHOLE_SEQ:
        results_df.to_csv(RESULTS_DIR / f'pklfiles/thdyinteg/fe_model_S288C_YAL002W-flanking1k_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt', index=False)
    else:
        results_df.to_csv(RESULTS_DIR / f'pklfiles/thdyinteg/fe_model_S288C_YAL002W-TSS_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt', index=False)



    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")