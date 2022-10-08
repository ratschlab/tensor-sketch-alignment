import re 
import csv
import time as mytime
import plotly.graph_objects as go
import json
from plotly.subplots import make_subplots
import argparse
import os
import subprocess
import numpy as np
from tqdm import tqdm
from pprint import pprint
from collections import defaultdict 
import edlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
#ASTARIX_PATH = "/home/alex/benchmark/datagen/astarix/release/astarix"
#VG_PATH = "/home/alex/benchmark/datagen/vg_source/vg/bin/vg"
#MINIGRAPH_PATH = "/home/alex/benchmark/datagen/minigraph/minigraph"
#DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"

NUM_THREADS = 8
METAGRAPH_PATH = "/cluster/apps/biomed/grlab/ameterez/metagraph/metagraph/build/metagraph"
VG_PATH = "/cluster/apps/biomed/grlab/ameterez/vg"
DATASET_DIR = None 
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]
THRESHOLDS = [0.05, 0.1]
NUM_QUERY_SEQS = 500

def rc(s):
    tab = str.maketrans("ACGT", "TGCA")
    return s.translate(tab)[::-1]

def vg_get_node_path_spelling(path, xg_file):
    pattern = ">\d+|\d+<"
    # >123>124>125
    nodes = re.findall(pattern, path.strip())
    full_spelling = ""
    reverse = False
    for node in nodes:
        import pdb
#        pdb.set_trace()
        if node[0] == '>':
            node = node[1:]
        elif node[-1] == '<':
            node = node[:-1]
            reverse=True
        command = f"{VG_PATH} find -x {xg_file} -n {node}"
        ret = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        command = f"{VG_PATH} view -j -"
        ret = subprocess.Popen(command.split(), stdin = ret.stdout, stdout = subprocess.PIPE)
        output = ret.communicate()[0]
#        print(output)
        spelling = json.loads(output)['node'][0]['sequence']
#        print(f"{node} -- {spelling.strip()}")
        full_spelling += rc(spelling.strip())
#    if reverse:
#        full_spelling = rc(full_spelling)
    return full_spelling
def read_seqs(path):
    lines = open(path, 'r').readlines()
    header_to_seq = {}
    n_lines = len(lines)
    for i in range(0, n_lines, 2):
        header = lines[i].strip()[1:]
        seq = lines[i+1].strip()
        header_to_seq[header] = seq

    return header_to_seq


def run_metagraph_default(mr):
    print("Running Default...")
    recalls = defaultdict(list) 
    times = []
    memory = []
    ref_seqs = read_seqs(os.path.join(DATASET_DIR, f'reference_{mr}.fa'))
    query_seqs = read_seqs(os.path.join(DATASET_DIR, f'mutated_{mr}.fa'))
    scores = np.ones(NUM_QUERY_SEQS, dtype=float) 
    command = f"/usr/bin/time -f %M {METAGRAPH_PATH} align " \
                  f"--seeder default " \
                  f"--parallel {NUM_THREADS} " \
                  f"--batch-size 1000 " \
                  f"--align-min-seed-length 15 " \
                  f"--align-xdrop 15 " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  f"{DATASET_DIR}/mutated_{mr}.fa "
    print(command)
    result_ = subprocess.run(command.split(), capture_output=True, text=True)
    memory.append(int(result_.stderr.strip().split("\n")[-1])/1024/1024)
    alignments = result_.stdout.split('\n')[:-2]
    time = json.loads(result_.stdout.split('\n')[-2])['time']
    avg_length = 0.0
    for alignment in alignments:
        alignment = alignment.split('\t')
        seq_idx = alignment[0]
        seq_idx_num = int(seq_idx[1:])
        matched_seq = alignment[3]
        ref_seq = ref_seqs[seq_idx]
        avg_length += len(matched_seq)
        edit_distance = edlib.align(ref_seq, matched_seq, task='path')['editDistance']
        scores[seq_idx_num] = edit_distance / len(ref_seq)
    times.append(time/NUM_QUERY_SEQS)
    return np.array(scores, dtype=float), np.array(times, dtype=float), np.array(memory, dtype=float)

def run_metagraph_sketching(mr, suffix):
    print("Running Sketching...")
    recalls = defaultdict(list) 
    times = []
    memory = []
    index_path = f"{DATASET_DIR}/sketch_mutation_{mr}_{suffix}.faiss"
    load_index = 1 if os.path.exists(index_path) else 0

    ref_seqs = read_seqs(os.path.join(DATASET_DIR, f'reference_{mr}.fa'))
    query_seqs = read_seqs(os.path.join(DATASET_DIR, f'mutated_{mr}.fa'))
    scores = np.ones(NUM_QUERY_SEQS, dtype=float) 
    command = f"/usr/bin/time -f %M {METAGRAPH_PATH} align " \
                  f"--seeder sketch " \
                  f"--embed-dim 14 " \
                  f"--n-times-sketch 1 " \
                  f"--parallel {NUM_THREADS} " \
                  f"--batch-size 1000 " \
                  f"--num-neighbours 10 " \
                  f"--align-end-bonus 0 " \
                  f"--index-path {index_path} " \
                  f"--load-index {load_index} " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  f"{DATASET_DIR}/mutated_{mr}.fa "
    #command = f"/usr/bin/time -f %M {METAGRAPH_PATH} align " \
    #              f"--seeder sketch " \
    #              f"--embed-dim 4 " \
    #              f"--n-times-sketch 1 " \
    #              f"--parallel {NUM_THREADS} " \
    #              f"--batch-size 1000 " \
    #              f"--num-neighbours 10 " \
    #              f"--align-end-bonus 0 " \
    #              f"-i {DATASET_DIR}/sequence_{K}.dbg " \
    #              f"{DATASET_DIR}/mutated_{mr}.fa "
    print(command)
    full_start_time = mytime.time()
    result_ = subprocess.run(command.split(), capture_output=True, text=True)
    full_end_time = mytime.time()
    memory.append(int(result_.stderr.strip().split("\n")[-1])/1024/1024)
    alignments = result_.stdout.split('\n')[:-2]
    time = json.loads(result_.stdout.split('\n')[-2])['time']
    for alignment in alignments:
        alignment = alignment.split('\t')
        seq_idx = alignment[0]
        seq_idx_num = int(seq_idx[1:])
        matched_seq = alignment[3]
        ref_seq = ref_seqs[seq_idx]
        query_seq = query_seqs[seq_idx]
        edit_distance = edlib.align(ref_seq, matched_seq, task='path')['editDistance']
        scores[seq_idx_num] = edit_distance / len(ref_seq)
    times.append(time/NUM_QUERY_SEQS)

    #with open(f"{DATASET_DIR}/time_output_{mr}.sketch", 'w') as timefile:
    #    timefile.write(str(full_end_time - full_start_time))

    return np.array(scores, dtype=float), np.array(times, dtype=float), np.array(memory, dtype=float)

def run_graphaligner(mr):
    print("Running graphaligner...")
    recalls = defaultdict(list) 
    times = []
    memory = []
    ref_seqs = read_seqs(os.path.join(DATASET_DIR, f'reference_{mr}.fa'))
    scores = np.ones(NUM_QUERY_SEQS, dtype=float) 
    #command = f"/usr/bin/time -f %M GraphAligner --precise-clipping 0.501 -t {NUM_THREADS} -g {DATASET_DIR}/sequence_{K}_blunted.gfa -f {DATASET_DIR}/mutated_{mr}.fa -a {DATASET_DIR}/aln_{mr}.gaf -x vg" #-x dbg"
    command = f"/usr/bin/time -f %M GraphAligner -t {NUM_THREADS} -g {DATASET_DIR}/sequence_{K}.gfa -f {DATASET_DIR}/mutated_{mr}.fa -a {DATASET_DIR}/aln_{mr}.gaf -x dbg"
    print(command) 
    time_start = mytime.time()
    result_ = subprocess.run(command.split(), capture_output=True, text=True)
    time_end = mytime.time()

    memory.append(int(result_.stderr.strip().split("\n")[-1])/1024/1024)
    command = f"python gfa_to_gaf.py {DATASET_DIR}/sequence_{K}.gfa {DATASET_DIR}/mutated_{mr}.fa {DATASET_DIR}/aln_{mr}.gaf {K}"
    result = subprocess.run(command.split(), capture_output=True, text=True)
    alignments = result.stdout.strip().split('>')[1:]
    for alignment in alignments:
        header, extra_path, score, cigar, path = alignment.strip().split(':')
        if path == 'N':
            continue
        seq_idx = header
        seq_idx_num = int(seq_idx[1:])
        ref_seq = ref_seqs[seq_idx]
        edit_distance = edlib.align(ref_seq, path, task='path')['editDistance'] 
        scores[seq_idx_num] = edit_distance / len(ref_seq)
    times.append((time_end-time_start) / NUM_QUERY_SEQS)
    return np.array(scores, dtype=float), np.array(times, dtype=float), np.array(memory, dtype=float)

def run_vg_mpmap(mr):
    print("Running vg mpmap...")
    recalls = defaultdict(list) 
    times = []
    memory = []

    ref_seqs = read_seqs(os.path.join(DATASET_DIR, f'reference_{mr}.fa'))
    query_seqs = read_seqs(os.path.join(DATASET_DIR, f'mutated_{mr}.fa'))
    
    scores = np.ones(NUM_QUERY_SEQS, dtype=float)
   
    f = open(f'{DATASET_DIR}/out_{mr}_mpmap.gaf','w')
    command = f"/usr/bin/time -f %M {VG_PATH} mpmap -t {NUM_THREADS} -z 1 -o 1 -n DNA -F GAF -x {DATASET_DIR}/sequence_map.xg -g {DATASET_DIR}/sequence_map.gcsa -f {DATASET_DIR}/mutated_{mr}.fa"
    print(command)
    g = open(f'{DATASET_DIR}/tmp_mpmap_{mr}', 'w')
    time_start = mytime.time()
    subprocess.run(command.split(), stderr=g, text=True, stdout=f)
    time_end = mytime.time()
    f.close()
    g.close()
    g = open(f'{DATASET_DIR}/tmp_mpmap_{mr}', 'r')
    result_ = g.read()
    g.close()
    memory.append(int(result_.strip().split("\n")[-1])/1024/1024)
    alignments = open(f"{DATASET_DIR}/out_{mr}_mpmap.gaf", "r").readlines()
    for alignment in alignments:
        parts = alignment.split("\t")
        if len(parts) < 13:
            continue
        start = int(parts[7])
        end = int(parts[8]) + 1
        path_spelling = vg_get_node_path_spelling(parts[5], f"{DATASET_DIR}/sequence_map.xg")[start:end]
        seq_idx = parts[0]
        seq_idx_num = int(seq_idx[1:])
        ref_seq = ref_seqs[seq_idx]
        query_seq = query_seqs[seq_idx]
        edit_distance = edlib.align(ref_seq, path_spelling, task='path')['editDistance']
        scores[seq_idx_num] = edit_distance / len(ref_seq)
    times.append((time_end-time_start)/NUM_QUERY_SEQS)
    return np.array(scores, dtype=float), np.array(times, dtype=float), np.array(memory, dtype=float)

def run_vg_map(mr):
    print("Running vg map...")
    recalls = defaultdict(list) 
    times = []
    memory = []
    ref_seqs = read_seqs(os.path.join(DATASET_DIR, f'reference_{mr}.fa'))
    query_seqs = read_seqs(os.path.join(DATASET_DIR, f'mutated_{mr}.fa'))
    scores = np.ones(NUM_QUERY_SEQS, dtype=float)
    f = open(f'{DATASET_DIR}/out_{mr}_map.gaf','w')
    command = f"/usr/bin/time -f %M {VG_PATH} map -t {NUM_THREADS} -z 1 -o 1 -x {DATASET_DIR}/sequence_map.xg -g {DATASET_DIR}/sequence_map.gcsa -f {DATASET_DIR}/mutated_{mr}.fa --gaf"
    print(command)
    g = open(f"{DATASET_DIR}/tmp_map_{mr}", 'w')
    time_start = mytime.time()
    subprocess.run(command.split(), stderr=g, text=True, stdout=f)
    time_end = mytime.time()
    f.close()
    g.close()
    g = open(f"{DATASET_DIR}/tmp_map_{mr}", 'r')
    result_ = g.read()
    g.close()
    memory.append(int(result_.strip().split("\n")[-1])/1024/1024)
    alignments = open(f"{DATASET_DIR}/out_{mr}_map.gaf", "r").readlines()
    for alignment in alignments:
        parts = alignment.split("\t")
        if len(parts) < 13:
            continue
        start = int(parts[7])
        end = int(parts[8]) + 1
        path_spelling = vg_get_node_path_spelling(parts[5], f"{DATASET_DIR}/sequence_map.xg")[start:end]
        seq_idx = parts[0]
        seq_idx_num = int(seq_idx[1:])
        ref_seq = ref_seqs[seq_idx]
        query_seq = query_seqs[seq_idx]
        edit_distance = edlib.align(ref_seq, path_spelling, task='path')['editDistance']
        scores[seq_idx_num] = edit_distance / len(ref_seq)
    times.append((time_end-time_start)/NUM_QUERY_SEQS)
    return np.array(scores, dtype=float), np.array(times, dtype=float), np.array(memory, dtype=float)


import pickle
from matplotlib import lines
sns.set_theme()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", type=str, required=True)
    parser.add_argument("--method", type=str, required=True)
    parser.add_argument("--mutation", type=int, required=True)
    parser.add_argument("--suffix", type=str, required=True)
    args = parser.parse_args()

    assert args.method in ['mpmap', 'map', 'default', 'sketch', 'graphaligner']

    DATASET_DIR = args.dataset_dir
    suffix = args.suffix
    #assert args.mutation in MUTATIONS

    method_scores = None
    method_times = None
    method_memory = None
    
    method_pkl = os.path.join(DATASET_DIR, f"{args.method}_{args.mutation}_{suffix}.pkl")
    if args.method == 'mpmap':
        method_scores, method_times, method_memory = run_vg_mpmap(args.mutation) 
    elif args.method == 'default':
        method_scores, method_times, method_memory = run_metagraph_default(args.mutation) 
    elif args.method == 'sketch':
        method_scores, method_times, method_memory = run_metagraph_sketching(args.mutation, suffix) 
    elif args.method == 'map':
        method_scores, method_times, method_memory = run_vg_map(args.mutation) 
    elif args.method == 'graphaligner':
        method_scores, method_times, method_memory = run_graphaligner(args.mutation) 
    # Save data 
    pickle.dump([method_scores, method_times, method_memory], open(method_pkl, 'wb')) 
