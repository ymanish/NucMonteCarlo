import os 
import numpy as np
from typing import Tuple
from py_analysis.config.gen_var import DATA_DIR, TI_PARAMS

def generate_parameter_file(mu_range:Tuple[float,float], 
                            sequences_file:str, 
                            parameter_file:str, 
                            n_mu:int=20)-> None:
    
    mu_values = np.round(np.linspace(mu_range[0], mu_range[1], n_mu), 2)

    with open(sequences_file, 'r') as f:
        entries = [line.strip().split('\t') for line in f.readlines()]
        positions = [int(entry[0]) for entry in entries]
        sequences = [entry[1] for entry in entries]

    if len(sequences) < 2:
        raise ValueError("The sequences file must contain at least two sequences.")

    ref_seq = sequences[0]
    ref_pos = positions[0]
    target_entries = list(zip(positions[1:], sequences[1:]))

    
    # first_sequence = "TTCCAAGGCCGATGAATTCGACTCTTTCCCAGCTGCCTCTGCTGCCGCTGCCGAAGAAGAAGAAGATGACGATGTCGATTTATTCGGTTCCGACGATGAAGAAGCTGACGCTGAAGCTGAAAAGTTGAAGGCTGAAAGAATTGCCGC"
    # other_sequences = "TATTAATGACACTCAAAAGACTTTCCTAGAATTTAGATCGTATACCCAATTAAGTGAAAAACTGGCATCTAGTTCTTCATATACGGCACCTCCCCTGAACGAAGATGGTCCTAAAGGGGTAGCTTCTGCAGTGTCACAAGGCTCCGA"
    # index = 1
    # with open(parameter_file, 'w') as f:
    #     for j, mu in enumerate(mu_values, start=1):
    #         f.write(f"{index}\t{mu:.2f}\t{first_sequence}\t{other_sequences}\n")
    #         index += 1
    
    with open(parameter_file, 'w') as f:
        # f.write("index\tmu\tseq_pos_start\tseq_ref\tseq_target\n")  # Header
        index = 1
        for pos, seq in target_entries:
            for mu in mu_values:
                f.write(f"{index}\t{mu:.2f}\t{pos}\t{ref_seq}\t{seq}\n")
                index += 1

    print(f"Parameter file generated at {parameter_file} with {len(sequences)-1} sequences.")


if __name__ == "__main__":
  

    # parameter_file = os.path.expanduser('~/pol/Projects/Codebase/NucleosomeMMC/State/Longer_Seq_SingleRepfe_PMF.txt')
    SEQ_LENGTH = TI_PARAMS['SEQ_LENGTH']
    SLIDE_INTERVAL = TI_PARAMS['SLIDE_INTERVAL']
    MU_RANGE = TI_PARAMS['MU_RANGE']
    MU_VALUES = TI_PARAMS['MU_VALUES']

    if TI_PARAMS['WHOLE_SEQ']:
        seq_list_file = DATA_DIR / f"S288C_YAL002W-flanking1k_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}.txt"
        parameter_file_path = DATA_DIR / f"S288C_YAL002W-flanking1k_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt"
    
    else:
        seq_list_file = DATA_DIR / f"S288C_YAL002W-TSS_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}.txt"
        parameter_file_path = DATA_DIR / f"S288C_YAL002W-TSS_win{SEQ_LENGTH}_sl{SLIDE_INTERVAL}_mu{MU_RANGE[0]}-{MU_RANGE[1]}_n{MU_VALUES}.txt"


    # Generate the parameter file
    generate_parameter_file(mu_range=MU_RANGE, sequences_file=seq_list_file, parameter_file=parameter_file_path, n_mu=MU_VALUES)







