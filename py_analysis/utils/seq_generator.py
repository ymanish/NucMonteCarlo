import random
from py_analysis.config.io_path import DATA_DIR
from typing import Dict
def generate_random_sequences(num_sequences: int, sequence_length: int, seed: int) -> dict:
 
    rng = random.Random(seed)
    bases = ["A", "T", "C", "G"]
    seq_dict = {}
    for i in range(num_sequences):
        seq = "".join(rng.choice(bases) for _ in range(sequence_length))
        seq_dict[f"seq_{i}"] = seq
    return seq_dict


def generate_gc_content_sequences(num_sequences: int, sequence_length: int, seed: int) -> Dict[str, str]:

    rng = random.Random(seed)
    seq_dict = {}

    # Calculate GC content step size (from 0% to 100%)
    if num_sequences == 1:
        gc_contents = [0.5]  
    else:
        gc_contents = [i / (num_sequences - 1) for i in range(num_sequences)]  # Evenly spaced from 0 to 1
    # print (gc_contents)
    for i, gc_content in enumerate(gc_contents):
        sequence = []
        for _ in range(sequence_length):
            # Probability of G or C is gc_content, A or T is (1 - gc_content)
            if rng.random() < gc_content:
                base = rng.choice(["G", "C"])
            else:
                base = rng.choice(["A", "T"])
            sequence.append(base)
        seq_dict[f"seq_{i}"] = "".join(sequence)

    return seq_dict



if __name__ == "__main__":
    num_sequences = 10000
    sequence_length = 147
    seed = 42
    output_file_path = DATA_DIR / f"gc_range_rand_seqs_{num_sequences//1000}k_seed_{seed}.txt"


    # random_sequences = generate_random_sequences(num_sequences, sequence_length, seed)
    # for key, value in random_sequences.items():
    #     print(f"{key}: {value}")

    gc_random_sequences = generate_gc_content_sequences(num_sequences, sequence_length, seed)

    with output_file_path.open("w") as f:
        for key, sequence in gc_random_sequences.items():
            f.write(f"{key}\t{sequence}\n")


    