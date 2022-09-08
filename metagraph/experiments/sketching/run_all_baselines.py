import re 
import csv
import time as mytime
import plotly.graph_objects as go
import multiprocessing
import json
from plotly.subplots import make_subplots
import argparse
import parasail
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

METAGRAPH_PATH = "/home/alex/metagraph/metagraph/build/metagraph"
ASTARIX_PATH = "/home/alex/benchmark/datagen/astarix/release/astarix"
VG_PATH = "/home/alex/benchmark/datagen/vg_source/vg/bin/vg"
MINIGRAPH_PATH = "/home/alex/benchmark/datagen/minigraph/minigraph"
DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
K = 80
MUTATIONS = [0]#, 5, 10, 15, 20, 25]
THRESHOLDS = [0.05, 0.1]
NUM_QUERY_SEQS = 500

def vg_get_node_path_spelling(path, xg_file):
    pattern = ">\d+|\d+<"
    # >123>124>125
    nodes = re.findall(pattern, path.strip())
    full_spelling = "" 
    for node in nodes:
        if node[0] == '>':
            node = node[1:]
        elif node[-1] == '<':
            node = node[:-1]
        command = f"{VG_PATH} find -x {xg_file} -n {node}"
        ret = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        command = f"{VG_PATH} view -j -"
        ret = subprocess.Popen(command.split(), stdin = ret.stdout, stdout = subprocess.PIPE)
        output = ret.communicate()[0]
        spelling = json.loads(output)['node'][0]['sequence']

        full_spelling += spelling.strip()
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

def init_scores(seqs):
    n_seqs = len(seqs)
    scores = []
    for key, val in seqs.items():
        scores.append(len(val))
    return np.array(scores, dtype=float)

def run_metagraph_default():
    print("Running Default...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        query_seqs = read_seqs(DATASET_DIR + f'mutated_{mr}.fa')
    
        scores = np.ones(NUM_QUERY_SEQS, dtype=float) 

        # command = f"{METAGRAPH_PATH} align " \
        #               f"--seeder default " \
        #               f"--parallel 23 " \
        #               f"--batch-size 1000 " \
        #               f"-i {DATASET_DIR}/sequence_{K}.dbg " \
        #               f"{DATASET_DIR}/mutated_{mr}.fa "
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder default " \
                      f"--parallel 23 " \
                      f"--batch-size 1000 " \
                      f"--align-xdrop 5 " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        print(command)
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        alignments = result_.stdout.split('\n')[:-2]
        time = json.loads(result_.stdout.split('\n')[-2])['time']
        avg_length = 0.0
        for alignment in alignments:
            print(alignment)
            import pdb
            # pdb.set_trace()
            alignment = alignment.split('\t')
            seq_idx = alignment[0]
            seq_idx_num = int(seq_idx[1:])
            matched_seq = alignment[3]
            ref_seq = ref_seqs[seq_idx]
            avg_length += len(matched_seq)
            edit_distance = edlib.align(ref_seq, matched_seq, task='path')['editDistance']
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append(time/NUM_QUERY_SEQS)
        print(times)
    return recalls, np.array(times, dtype=float)

def run_metagraph_sketching():
    print("Running Sketching...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        query_seqs = read_seqs(DATASET_DIR + f'mutated_{mr}.fa')
    
        scores = np.ones(NUM_QUERY_SEQS, dtype=float) 

        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder sketch " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--parallel 23 " \
                      f"--batch-size 1000 " \
                      f"--num-neighbours 10 " \
                      f"--align-end-bonus 0 " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        # command = f"{METAGRAPH_PATH} align " \
        #               f"--seeder sketch " \
        #               f"--embed-dim 14 " \
        #               f"--n-times-sketch 1 " \
        #               f"--parallel 23 " \
        #               f"--batch-size 1000 " \
        #               f"--num-neighbours 10 " \
        #               f"--align-edit-distance " \
        #               f"--align-end-bonus 0 " \
        #               f"-i {DATASET_DIR}/sequence_{K}.dbg " \
        #               f"{DATASET_DIR}/mutated_{mr}.fa "
        print(command)
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
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

        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append(time/NUM_QUERY_SEQS)
        print(times)
    return recalls, np.array(times, dtype=float)

def run_graphaligner():
    print("Running graphaligner...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float) 
       
        # command = f"GraphAligner --precise-clipping 0.501 -g data/sequence_{K}.gfa -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -x dbg" #-x dbg"
        command = f"GraphAligner --precise-clipping 0.501 -g data/sequence_{K}_blunted.gfa -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -x vg" #-x dbg"
        print(command) 
        time_start = mytime.time()
        subprocess.run(command.split())
        time_end = mytime.time()

        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            seq_idx = alignment[0]
            seq_idx_num = int(seq_idx[1:])
            ref_seq = ref_seqs[seq_idx]
            edit_distance = float(alignment[12][5:])
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((time_end-time_start) / NUM_QUERY_SEQS)
        print(times)
    return recalls, np.array(times, dtype=float)

def run_vg_mpmap():
    print("Running vg mpmap...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        query_seqs = read_seqs(DATASET_DIR + f'mutated_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float)
        f = open(f'data/out_{mr}_mpmap.gaf','w')
        command = f"{VG_PATH} mpmap -z 1 -o 1 -n DNA -F GAF -x data/sequence_map.xg -g data/sequence_map.gcsa -f data/mutated_{mr}.fa"
        print(command)
        time_start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        time_end = mytime.time()
        f.close()

        alignments = open(f"data/out_{mr}_mpmap.gaf", "r").readlines()
        for alignment in alignments:
            import pdb
            # pdb.set_trace()
            parts = alignment.split("\t")
            if len(parts) < 13:
                continue
            start = int(parts[7])
            end = int(parts[8])
            path_spelling = vg_get_node_path_spelling(parts[5], "data/sequence_map.xg")[start:end]
            seq_idx = parts[0]
            seq_idx_num = int(seq_idx[1:])
            ref_seq = ref_seqs[seq_idx]
            query_seq = query_seqs[seq_idx]
            edit_distance = edlib.align(ref_seq, path_spelling, task='path')['editDistance']
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((time_end-time_start)/NUM_QUERY_SEQS)
        print(times)
    return recalls, np.array(times, dtype=float)

def run_vg_map():
    print("Running vg map...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        query_seqs = read_seqs(DATASET_DIR + f'mutated_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float)
        f = open(f'data/out_{mr}_map.gaf','w')
        command = f"{VG_PATH} map -z 1 -o 1 -x data/sequence_map.xg -g data/sequence_map.gcsa -f data/mutated_{mr}.fa --gaf"
        print(command)
        time_start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        time_end = mytime.time()
        f.close()

        alignments = open(f"data/out_{mr}_map.gaf", "r").readlines()
        for alignment in alignments:
            parts = alignment.split("\t")
            if len(parts) < 13:
                continue
            start = int(parts[7])
            end = int(parts[8])
            path_spelling = vg_get_node_path_spelling(parts[5], "data/sequence_map.xg")[start:end]
            seq_idx = parts[0]
            seq_idx_num = int(seq_idx[1:])
            ref_seq = ref_seqs[seq_idx]
            query_seq = query_seqs[seq_idx]
            edit_distance = edlib.align(ref_seq, path_spelling, task='path')['editDistance']
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((time_end-time_start)/NUM_QUERY_SEQS)
        print(times)
    return recalls, np.array(times, dtype=float)


import pickle
from matplotlib import lines
sns.set_theme()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--regen", action="store_true")
    args = parser.parse_args()
    if args.regen:
        print("REGENERATING SEQUENCES")
        # vg_mpmap_recalls, vg_mpmap_times = run_vg_mpmap() 
        # vg_map_recalls, vg_map_times = run_vg_map() 
        # sketch_recalls, sketch_times = run_metagraph_sketching()
        default_recalls, default_times = run_metagraph_default()
        # graphaligner_recalls, graphaligner_times = run_graphaligner()

        # Save data 
        # pickle.dump([vg_mpmap_recalls, vg_mpmap_times], open("mpmap.pkl", 'wb')) 
        # pickle.dump([vg_map_recalls, vg_map_times], open("map.pkl", 'wb')) 
        # pickle.dump([sketch_recalls, sketch_times], open("sketch.pkl", 'wb')) 
        # pickle.dump([default_recalls, default_times], open("default.pkl", 'wb')) 
        # pickle.dump([graphaligner_recalls, graphaligner_times], open("graphaligner.pkl", 'wb')) 
    else:
        print("USING OLD SEQUENCES")
    
    linestyles = list(lines.lineStyles.keys())
    vg_mpmap_recalls, vg_mpmap_times = pickle.load(open("mpmap.pkl", "rb"))
    vg_map_recalls, vg_map_times = pickle.load(open("map.pkl", "rb"))
    sketch_recalls, sketch_times = pickle.load(open("sketch.pkl", "rb"))
    default_recalls, default_times = pickle.load(open("default.pkl", "rb"))
    graphaligner_recalls, graphaligner_times = pickle.load(open("graphaligner.pkl", "rb"))
    
    plt.figure(figsize=(20,10))
    for T_idx, T in enumerate(THRESHOLDS):
        plt.plot(MUTATIONS, vg_mpmap_recalls[T], color='r', label=f"vg mpmap @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, vg_map_recalls[T], color='g', label=f"vg map @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, sketch_recalls[T], color='b', label=f"Sketch @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, default_recalls[T], color='y', label=f"Default @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, graphaligner_recalls[T], color='m', label=f"GraphAligner @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
    plt.ylabel("Recall")
    plt.xlabel("Mutation rate")
    plt.title("Baselines")
    plt.legend()
    plt.savefig("mutation_baselines.pdf")
    plt.close()

    plt.figure(figsize=(20,10))
    plt.plot(MUTATIONS, vg_mpmap_times, color='r', label=f"vg mpmap", marker='o', markersize=6)
    plt.plot(MUTATIONS, vg_map_times, color='g', label=f"vg map", marker='o', markersize=6)
    plt.plot(MUTATIONS, sketch_times, color='b', label=f"Sketch", marker='o', markersize=6)
    plt.plot(MUTATIONS, default_times, color='y', label=f"Default", marker='o', markersize=6)
    plt.plot(MUTATIONS, graphaligner_times, color='m', label=f"GraphAligner", marker='o', markersize=6)
    plt.yscale("log")
    plt.ylabel("Time")
    plt.xlabel("Mutation rate")
    plt.title("Baselines")
    plt.legend()
    plt.savefig("time_baselines.pdf")
    plt.close()
    print("Done")
