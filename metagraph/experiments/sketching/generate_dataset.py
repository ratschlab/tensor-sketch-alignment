import os
import numpy as np
import subprocess
import random
import shutil
import argparse

def get_random_str(main_str, substr_len):
    idx = random.randrange(0, len(main_str) - substr_len + 1)
    return main_str[idx:(idx+substr_len)]


DATASET_DIR = "./data"

get_blunted_path = "/home/alex/benchmark/datagen/GetBlunted/build/get_blunted"
vg_path = "/home/alex/benchmark/datagen/vg"
print(vg_path)
print(get_blunted_path)
# Clean up
shutil.rmtree(DATASET_DIR)
os.mkdir(DATASET_DIR)

INPUT_SEQ = "sequence.fa"
ALPHABET = ['A', 'C', 'T', 'G']
if __name__ == '__main__':
    assert os.path.exists(vg_path)
    assert os.path.exists(get_blunted_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--metagraph-path", type=str, required=True, help="Path to metagraph executable")
    parser.add_argument("--graph-seq-len", type=int, required=True, help="Length of the seq that the graph is generated from")
    parser.add_argument("--max-k", type=int, required=True, help="Maximum k-mer size (up to 85)")
    args = parser.parse_args()

    METAGRAPH_PATH = args.metagraph_path
    GRAPH_SEQ_LEN = args.graph_seq_len
    MAX_K = args.max_k

    assert os.path.exists(DATASET_DIR), "Please create dataset directory"
    assert MAX_K < GRAPH_SEQ_LEN, "Choose a smaller K"
    
    # Generate a random sequence
    graph_seq = ''.join(np.random.choice(ALPHABET, GRAPH_SEQ_LEN))
    header = f'>Base sequence'
    out_graph_seq = '\n'.join([header, graph_seq])
    # print("[LOG] Generated graph sequence:")
    # print(out_graph_seq)
    graph_seq_path = os.path.join(DATASET_DIR, INPUT_SEQ)
    with open(graph_seq_path, 'w') as f:
        f.write(out_graph_seq)

    # Generate graph
    for K in range(20, MAX_K, 10):
        dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}')
        blunted_dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}_blunted')
        build_command = f"{METAGRAPH_PATH} build -k {K} -o {dbg_output}.dbg {graph_seq_path}"
        assemble_command = f"{METAGRAPH_PATH} assemble --to-gfa --compacted --unitigs -o {dbg_output}.gfa {dbg_output}.dbg"
        blunt_command = f"{get_blunted_path} --input_gfa {dbg_output}.gfa"
        vg_command = f"{vg_path} -g blunted_{dbg_output}.gfa -v"
        subprocess.run(build_command.split())
        subprocess.run(assemble_command.split())
        open(f"{blunted_dbg_output}.gfa", 'w').write(str(subprocess.run(blunt_command.split(), capture_output=True).stdout))
        open(f"{blunted_dbg_output}.vg", 'w').write(str(subprocess.run(vg_command.split(), capture_output=True).stdout))
        print(f"[LOG] Saved .dbg file from generated sequence - {K}")
