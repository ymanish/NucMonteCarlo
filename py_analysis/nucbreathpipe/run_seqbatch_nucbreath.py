
##This script is used to run the nucleosome breathing energies in parallel for mutliple number of sequences.
##    Takes the KRESCFACTOR as input from the command line.
#     STYLE = "ph_index" # "b_index" or "open_sites" or "ph_index"
#     LENGTH = 28 # 28 for ph_index, 14 for b_index
#     nucmethod = "hybrid"
#     FOR_DNA_HISTONE = True  # if you use ph_index, set this to True because it saves the data in the dna_histone folder, otherwise it saves in the nucbreathfe folder
#     FLIP_SEQUENCE = False
#     HARD_CONS = False
#     # seq_dict needs to be given for the sequences the breathing energies are calculated for.

from nucbreath_calculator import BreathEnergies_per_sequence

import sys
import os
from typing import List
import concurrent.futures
from tqdm import tqdm
import pandas as pd
from py_analysis.config.breath_var import Parameters, seq_dict, MAX_WORKERS
from breathstates_calculator import total_states_index, total_open_states
from process_breathresults import process_results, save_results
from py_analysis.config.custom_types import NuclBreathingResult
from py_analysis.config.io_path import DATA_DIR
import argparse
from pathlib import Path


def slrum_job_fn(batch_save_dir_prefix:str):
    # Configure argument parser
    parser = argparse.ArgumentParser(
        description="Process nucleosome breathing states for a sequence batch for the SLURM cluster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--kresc", type=float, required=True,
                       help="KRESCFACTOR value")
    parser.add_argument("--start", type=int, default=0,
                       help="Start index (inclusive)")
    parser.add_argument("--end", type=int, required=True,
                       help="End index (exclusive)")
    parser.add_argument("--batch", type=int, required=True,
                       help="Batch/SLURM array task ID")
    parser.add_argument("--seq-file", type=Path, required=True,
                       help="Path to sequence file, the file should be in the format: id\tsequence")

    args = parser.parse_args()

    # Validate indices
    if args.start < 0:
        raise ValueError(f"Start index cannot be negative: {args.start}")
    if args.end <= args.start:
        raise ValueError(f"End index ({args.end}) must be > start ({args.start})")

    print(f"\nBatch {args.batch} parameters:")
    print(f"KRESCFACTOR: {args.kresc}")
    print(f"Processing sequences: {args.start}-{args.end-1}")

    # Load sequences
    try:
        seq_df = pd.read_csv(args.seq_file, sep="\t", 
                            header=None, names=["id", "sequence"])
    except FileNotFoundError:
        raise SystemExit(f"ERROR: Sequence file not found: {args.seq_file}")

    if args.end > len(seq_df):
        raise ValueError(f"End index ({args.end}) exceeds total sequences ({len(seq_df)})")

    # Process batch
    batch_df = seq_df.iloc[args.start:args.end]
    seq_dict = dict(zip(batch_df['id'], batch_df['sequence']))
    
    save_dir = Path(f"{batch_save_dir_prefix}_{args.batch}")
        
    return seq_dict, save_dir, args.kresc

if __name__ == "__main__":
    
    import time 
    start = time.perf_counter()

    params = Parameters(KRESCFACTOR=float(sys.argv[1]))
    SAVE_DIR_STR = "nbfiles"

    ###### START HERE: FOR SlURM CLUSTER USE ##################
    ###########################################################
    ###########################################################
    
    # seq_dict, SAVE_DIR_STR, kfact = slrum_job_fn(batch_save_dir_prefix="gcseqs_breath_batch")
    # params = Parameters(KRESCFACTOR=float(kfact))
   

    ###### ENDS HERE: FOR SlURM CLUSTER USE ##################
    ###########################################################
    ###########################################################


    # Calculate binding states based on STYLE
    if params.STYLE == "open_sites":
        states = total_open_states()
    elif params.STYLE == "b_index":
        states = total_states_index(length=params.LENGTH)
    else:
        states = total_states_index(length=params.LENGTH)

    # Print parameters for logging
    print(f"PARAMETERS...............")
    print(f"nucmethod: {params.nucmethod}")
    print(f"Free DNA method: {params.freedna_method}")
    print(f"FOR DNA-HISTONE: {params.FOR_DNA_HISTONE}")
    print(f"Sequence Flipped: {params.FLIP_SEQUENCE}")
    print(f"Style: {params.STYLE}")
    print(f"Length: {params.LENGTH}")
    print(f"Factor: {params.KRESCFACTOR}")
    print(f"Hard constraints: {params.HARD_CONS}")
    print(f"Total states: {len(states)}")
    print(f"Total sequences: {len(seq_dict)}")

    results_all: List[NuclBreathingResult] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [
            executor.submit(
                BreathEnergies_per_sequence,
                key=s,
                seq=seq_dict[s],
                nucmethod=params.nucmethod,
                freedna_method=params.freedna_method,
                bind_sates=states,
                factor=params.KRESCFACTOR,
                hard=params.HARD_CONS,
                style=params.STYLE,
                flip=params.FLIP_SEQUENCE
            ) for s in seq_dict.keys()
        ]
        total = len(futures)
        for future in tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences"):
            energy_array = future.result()
            results_all.extend(energy_array)

    # Process and save results
    df_full = process_results(results_all, params)
    save_results(df_full, params, dir=SAVE_DIR_STR)

    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")