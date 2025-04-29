
import os
from py_analysis.config.gen_var import DATA_DIR, TI_PARAMS

def from_fasta_to_sliding_list(fasta_file:str, 
                                 output_file:str, 
                                 sequence_length:int=147, 
                                 slide_interval:int=1) -> None:
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))


    with open(output_file, 'w') as f:
        for i in range(0, len(sequence) - sequence_length + 1, slide_interval):
            running_sequence = sequence[i:i + sequence_length]
            f.write(f"{i}\t{running_sequence}\n") ## here i is the start index of the sequence

    return print(f"Sliding window sequences saved to {output_file}")

if __name__ == "__main__":

    SLIDE_INTERVAL = TI_PARAMS['SLIDE_INTERVAL']
    SEQUENCE_LENGTH = TI_PARAMS['SEQ_LENGTH']
    WHOLE_SEQ = TI_PARAMS['WHOLE_SEQ']


    if WHOLE_SEQ:

        #Uses the whole 5000 bp sequence basicially the whole gene with +1kbp upstream and +1kbp downstream)
        fasta_file_path = DATA_DIR / "S288C_YAL002W-flanking1k.fsa"
        output_file_path = DATA_DIR / f"S288C_YAL002W-flanking1k_win{SEQUENCE_LENGTH}_sl{SLIDE_INTERVAL}.txt"
    else:
        #Uses the shorter version of the gene, -500bp and +1000bp around the TSS)
        fasta_file_path = DATA_DIR / "S288C_YAL002W-TSS_region.fsa"
        output_file_path = DATA_DIR / f"S288C_YAL002W-TSS_win{SEQUENCE_LENGTH}_sl{SLIDE_INTERVAL}.txt"


 
    from_fasta_to_sliding_list(fasta_file_path, output_file_path, sequence_length=SEQUENCE_LENGTH, slide_interval=SLIDE_INTERVAL)




    ### changes i want to make:
    # 1. give the main fasta file as input and I get the # output_file = os.path.expanduser('~/pol/Projects/Codebase/NucleosomeMMC/State/S288C_YAL002W_selected_sequences.txt')
#or output_file = os.path.expanduser('~/pol/Projects/Codebase/NucleosomeMMC/State/YAL002W_selected_sequences.txt') (for shorter version) as direcly the output which I can feed to the gen_thermody_params.py
