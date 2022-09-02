import os
import numpy as np
import subprocess
import random
import shutil
import argparse

def get_random_str(main_str, substr_len):
    idx = random.randrange(0, len(main_str) - substr_len + 1)
    return main_str[idx:(idx+substr_len)]

def make_vg():
    command1 = f"{vg_path} construct -r data/sequence.fa"
    f = open("data/sequence.vg", 'w')
    subprocess.run(command1.split(), stdout=f)
    command2 = f"{vg_path} index -x data/sequence.xg -g data/sequence.gcsa -k 16 data/sequence.vg"
    subprocess.run(command2.split())

DATASET_DIR = "./data"
get_blunted_path = "/home/alex/benchmark/datagen/GetBlunted/build/get_blunted"
vg_path = "/home/alex/benchmark/datagen/vg"

# Clean up
# shutil.rmtree(DATASET_DIR)
# os.mkdir(DATASET_DIR)

INPUT_SEQ = "sequence.fa"
GENOME_SEQ = None 
ALPHABET = ['A', 'C', 'T', 'G']

def mutate(s, rate):
    mutated_string = ""
    i = 0
    while i < len(s):
        chance = np.random.random()
        if chance < rate:
            mutated_string += np.random.choice(['A', 'C', 'T', 'G'])
        else:
            mutated_string += s[i]
        i += 1
    return mutated_string

if __name__ == '__main__':
    assert os.path.exists(vg_path)
    assert os.path.exists(get_blunted_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--metagraph-path", type=str, required=True, help="Path to metagraph executable")
    parser.add_argument("--graph-seq-len", type=int, required=True, help="Length of the seq that the graph is generated from")
    parser.add_argument("--num-levels", type=int, required=True, help="Number of seqs")
    parser.add_argument("--mutation-rate", type=float, required=True, help="Mutation rate of the sequences")
    args = parser.parse_args()

    METAGRAPH_PATH = args.metagraph_path
    GRAPH_SEQ_LEN = args.graph_seq_len
    K = 80 
    graph_seq_path = os.path.join(DATASET_DIR, INPUT_SEQ)
    dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}')
    blunted_dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}_blunted')
    print(f"Levels: {args.num_levels + 1}")
    print(f"Mutation rate: {args.mutation_rate}")
    print(f"Kmer: {K}")
    print(f"Genome sequence: {GENOME_SEQ}")
    print(f"Graph seq path: {graph_seq_path}")

    assert os.path.exists(DATASET_DIR), "Please create dataset directory"
    
    if GENOME_SEQ is not None:
        graph_seq = open(os.path.join(DATASET_DIR, GENOME_SEQ), 'r').read()
    else:
        print("Generating random string because GENOME_SEQ is None")
        graph_seq = "".join(np.random.choice(['A', 'C', 'T', 'G'], args.graph_seq_len))
    seqs = [graph_seq]
    
    for i in range(args.num_levels):
        new_seqs = []
        for seq in seqs:
            new_s = mutate(seq, args.mutation_rate)
            new_seqs.append(new_s)
        seqs += new_seqs
    print(f"#Sequences: {len(seqs)}")
    seq_file = []
    for i in range(len(seqs)):
        header = f">Sequence{i}"
        seq_file += [header, seqs[i]]

    seq_output = '\n'.join(seq_file).strip()
    with open(os.path.join(DATASET_DIR, INPUT_SEQ), 'w') as f:
        f.write(seq_output)
    print("Done with the sequences")
    # make_vg()
   
    # Generate graph
    build_command = f"{METAGRAPH_PATH} build -k {K} --parallel 20 -o {dbg_output}.dbg {graph_seq_path}"
    print(build_command)
    subprocess.run(build_command.split())
    print(f"[LOG] Saved .dbg file from generated sequence - {K}")

    # Assemble
    assemble_command = f"{METAGRAPH_PATH} assemble --to-gfa --compacted --unitigs -o {dbg_output}.gfa {dbg_output}.dbg"
    subprocess.run(assemble_command.split())
   
    # Blunt graph
    blunt_command = f"{get_blunted_path} --input_gfa {dbg_output}.gfa"
    blunted_graph = subprocess.run(blunt_command.split(), capture_output=True).stdout.decode("utf-8")
    open(f"{blunted_dbg_output}.gfa", 'w').write(blunted_graph)
   
    # Build vg index
    vg_command = f"{vg_path} autoindex -g data/sequence_{K}_blunted.gfa -V 2 -w map --prefix data/sequence"
    subprocess.run(vg_command.split())
    
    print("Done")
