import os
import numpy as np
import subprocess
import random


def get_random_str(main_str, substr_len):
    idx = random.randrange(0, len(main_str) - substr_len + 1)
    return main_str[idx:(idx+substr_len)]


def mutate(s, alphabet, mutation_rate):
    new_s = []
    for i in range(len(s)):
        chance = random.random()
        if chance > mutation_rate: # no mutation
            new_s.append(s[i])
        elif chance <= mutation_rate: # mutate
            filtered_alphabet = [x for x in alphabet if x != s[i]]
            c_ = np.random.choice(filtered_alphabet)
            new_s.append(c_)
    return ''.join(new_s)

def mutate2(s, alphabet):
    return ''.join(np.random.choice(alphabet, len(s)))

DATASET_DIR = "./data"
INPUT_SEQ = "sequence.fa"
QUERY_SEQS = "query_strings.fa"
METAGRAPH_PATH = "/Users/alex/metagraph/metagraph/build/metagraph"
ALPHABET = ['A', 'C', 'T', 'G']

GRAPH_SEQ_LEN = 10
MUTATION_RATE = 0.0

NUM_QUERY_SEQS = 10
QUERY_STRING_LEN = 8
MAX_K = QUERY_STRING_LEN

if __name__ == '__main__':
    assert os.path.exists(DATASET_DIR), "Please create dataset directory"
    # assert K < GRAPH_SEQ_LEN, "Choose a smaller K"
    assert QUERY_STRING_LEN < GRAPH_SEQ_LEN, "Query string too long"

    # Generate a random sequence
    graph_seq = ''.join(np.random.choice(ALPHABET, GRAPH_SEQ_LEN))
    header = f'>Base sequence'
    out_graph_seq = '\n'.join([header, graph_seq])
    print("[LOG] Generated graph sequence:")
    print(out_graph_seq)
    graph_seq_path = os.path.join(DATASET_DIR, INPUT_SEQ)
    with open(graph_seq_path, 'w') as f:
        f.write(out_graph_seq)

    # Generate graph
    for K in range(2, MAX_K):
        build_command = f"{METAGRAPH_PATH} build -k {K} -o {os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}')} {graph_seq_path}"
        subprocess.run(build_command.split())
        print(f"[LOG] Saved .dbg file from generated sequence - {K}")

    # Generate query seqs
    query_strings = []
    for i in range(NUM_QUERY_SEQS):
        query_string = get_random_str(graph_seq, QUERY_STRING_LEN)
        query_string_mutated = mutate(query_string, alphabet=ALPHABET, mutation_rate=MUTATION_RATE)
        #query_string_mutated = mutate2(query_string, alphabet=ALPHABET)
        header = f">Q{i} - mutation rate {MUTATION_RATE}"
        query_strings.append('\n'.join([header, query_string_mutated]))
    query_strings_out = '\n'.join(query_strings)
    print(f"[LOG] Generated {NUM_QUERY_SEQS} query strings with mutation rate {MUTATION_RATE}")
    query_strings_path = os.path.join(DATASET_DIR, QUERY_SEQS)
    with open(query_strings_path, "w") as f:
        f.write(query_strings_out)
    print(f"[LOG] Saved query strings")

