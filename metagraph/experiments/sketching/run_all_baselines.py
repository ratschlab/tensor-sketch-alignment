import re 
import cigar
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
from cigar import Cigar
import numpy as np
from tqdm import tqdm
from pprint import pprint
from collections import defaultdict 
import edlib
import matplotlib.pyplot as plt
import seaborn as sns
METAGRAPH_PATH = "/home/alex/metagraph/metagraph/build/metagraph"
ASTARIX_PATH = "/home/alex/benchmark/datagen/astarix/release/astarix"
VG_PATH = "/home/alex/benchmark/datagen/vg_source/vg/bin/vg"
MINIGRAPH_PATH = "/home/alex/benchmark/datagen/minigraph/minigraph"
DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]
# MUTATIONS = [5]#, 10, 15, 20, 25]
# MUTATIONS = [0]
# THRESHOLDS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
THRESHOLDS = [0.05, 0.1, 0.2]
NUM_QUERY_SEQS = 300

def convert_to_cigar(s):
    # - deletion
    # + insertion
    # :3 matches
    # * substitution
    cigar_string = ""
    regex_expr = "-[AGCT]+|\+[AGCT]+|:[0-9]+|\*[ACGT]+"
    # regex_expr = "/(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/"
    matches = re.findall(regex_expr, s)
    for matchh in matches:
        op = matchh[0]
        arg = matchh[1:]
        if op == "-":
            cigar_string += f"{len(arg)}D"
        elif op == "+":
            cigar_string += f"{len(arg)}I"
        elif op == ":":
            cigar_string += f"{int(arg)}="
        elif op == "*":
            cigar_string += f"1X"
    return cigar_string

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

def ed_from_cigar(s):
    score = 0
    parts = list(cigar.Cigar(s).items())
    n_parts = len(parts)
    i = 0
    while i < n_parts:
        count, op = parts[i]
        if op == '=':
            score += count * 0 
        elif op == 'X':
            score += count * 1 
        elif op == 'I' or op == 'D':
            score += count * 1 
        i += 1
    return score
def score_from_cigar(s):
    score = 0
    parts = list(cigar.Cigar(s).items())
    n_parts = len(parts)
    i = 0
    while i < n_parts:
        count, op = parts[i]
        if op == '=':
            score += count * 0 
        elif op == 'X':
            score -= count * 1 
        elif op == 'I' or op == 'D':
            score -= (1 + (count-1) * 1)
        i += 1
    return score

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
                      f"--parallel 20 " \
                      f"--batch-size 1000 " \
                      f"--num-neighbours 10 " \
                      f"--align-edit-distance " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        print(command)
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        alignments = result_.stdout.split('\n')[:-2]
        time = json.loads(result_.stdout.split('\n')[-2])['time']
        for alignment in alignments:
            alignment = alignment.split('\t')
            
            seq_idx = alignment[0]
            seq_idx_num = int(seq_idx[1:])
            matched_seq = alignment[3]
            ref_seq = ref_seqs[seq_idx]

            edit_distance = edlib.align(ref_seq, matched_seq, task='path')['editDistance']

            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append(time/NUM_QUERY_SEQS)
    return recalls, np.array(times, dtype=float)

def run_graphaligner():
    print("Running graphaligner...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float) 
       
        command = f"GraphAligner -g data/sequence_{K}.gfa -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -x dbg" #-x dbg"
        command = f"GraphAligner -g data/sequence_{K}_blunted.gfa -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -x vg" #-x dbg"
        
        start = mytime.time()
        subprocess.run(command.split())
        end = mytime.time()

        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            print(alignment)
            seq_idx = alignment[0]
            seq_idx_num = int(seq_idx[1:])
            ref_seq = ref_seqs[seq_idx]
            edit_distance = float(alignment[12][5:])
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        print(scores)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((end-start) / NUM_QUERY_SEQS)
    return recalls, np.array(times, dtype=float)

def run_vg_mpmap():
    print("Running vg mpmap...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float)
        f = open(f'data/out_{mr}.gam','w')
        command = f"{VG_PATH} mpmap --match 1 --full-l-bonus 0 --mismatch 1 --gap-open 1 --gap-extend 1 -n DNA -F GAF -x data/sequence.xg -g data/sequence.gcsa -f data/mutated_{mr}.fa"
        start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        end = mytime.time()
        f.close()

        alignments = open(f"data/out_{mr}.gam", "r").readlines()
        for alignment in alignments:
            parts = alignment.split("\t")
            if len(parts) < 13:
                continue
            seq_idx = parts[0]
            seq_idx_num = int(seq_idx[1:])
            
            csz_string = parts[13][5:]
            cigar_str = convert_to_cigar(csz_string)
            
            ref_seq = ref_seqs[seq_idx]
            edit_distance = float(ed_from_cigar(cigar_str))
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((end-start)/NUM_QUERY_SEQS)
    return recalls, np.array(times)

def run_vg_map():
    print("Running vg map...")
    recalls = defaultdict(list) 
    times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = read_seqs(DATASET_DIR + f'reference_{mr}.fa')
        scores = np.ones(NUM_QUERY_SEQS, dtype=float)
        f = open(f'data/out_{mr}.gam','w')
        command = f"{VG_PATH} map --match 1 --full-l-bonus 0 --mismatch 1 --gap-open 1 --gap-extend 1 -x data/sequence.xg -g data/sequence.gcsa -f data/mutated_{mr}.fa --gaf"
        start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        end = mytime.time()
        f.close()

        alignments = open(f"data/out_{mr}.gam", "r").readlines()
        for alignment in alignments:
            parts = alignment.split("\t")
            if len(parts) < 13:
                continue
            seq_idx = parts[0]
            seq_idx_num = int(seq_idx[1:])
            
            csz_string = parts[13][5:]
            cigar_str = convert_to_cigar(csz_string)
            
            ref_seq = ref_seqs[seq_idx]
            edit_distance = float(ed_from_cigar(cigar_str))
            scores[seq_idx_num] = edit_distance / len(ref_seq)
        for threshold in THRESHOLDS:
            recalls[threshold].append(np.mean(scores <= threshold))
        times.append((end-start)/NUM_QUERY_SEQS)
    return recalls, np.array(times)


import pickle
from matplotlib import lines
sns.set_theme()
if __name__ == '__main__':
    linestyles = list(lines.lineStyles.keys())
    generate = False 
    if generate:
        vg_mpmap_recalls, vg_mpmap_times = run_vg_mpmap() 
        vg_map_recalls, vg_map_times = run_vg_map() 
        sketch_recalls, sketch_times = run_metagraph_sketching()
        graphaligner_recalls, graphaligner_times = run_graphaligner()

        # Save data 
        pickle.dump([vg_mpmap_recalls, vg_mpmap_times], open("mpmap.pkl", 'wb')) 
        pickle.dump([vg_map_recalls, vg_map_times], open("map.pkl", 'wb')) 
        pickle.dump([sketch_recalls, sketch_times], open("sketch.pkl", 'wb')) 
        pickle.dump([graphaligner_recalls, graphaligner_times], open("graphaligner.pkl", 'wb')) 

    
    vg_mpmap_recalls, vg_mpmap_times = pickle.load(open("mpmap.pkl", "rb"))
    vg_map_recalls, vg_map_times = pickle.load(open("map.pkl", "rb"))
    sketch_recalls, sketch_times = pickle.load(open("sketch.pkl", "rb"))
    graphaligner_recalls, graphaligner_times = pickle.load(open("graphaligner.pkl", "rb"))
    
    plt.figure(figsize=(20,10))
    for T_idx, T in enumerate(THRESHOLDS):
        plt.plot(MUTATIONS, vg_mpmap_recalls[T], color='r', label=f"vg mpmap @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, vg_map_recalls[T], color='g', label=f"vg map @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, sketch_recalls[T], color='b', label=f"Sketch @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(MUTATIONS, graphaligner_recalls[T], color='m', label=f"GraphAligner @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
    plt.ylabel("Recall")
    plt.xlabel("Mutation rates")
    plt.title("Baselines")
    plt.legend()
    plt.savefig("mutation_baselines.png")
    plt.close()

    plt.figure(figsize=(20,10))
    for T_idx, T in enumerate(THRESHOLDS):
        plt.plot(vg_mpmap_times, vg_mpmap_recalls[T], color='r', label=f"vg mpmap @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(vg_map_times, vg_map_recalls[T], color='g', label=f"vg map @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(sketch_times, sketch_recalls[T], color='b', label=f"Sketch @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
        plt.plot(graphaligner_times, graphaligner_recalls[T], color='m', label=f"GraphAligner @ {T}", linestyle=linestyles[T_idx], marker='o', markersize=6)
    plt.ylabel("Recall")
    plt.xlabel("Average time(s)")
    plt.title("Baselines")
    plt.legend()
    plt.savefig("time_baselines.png")
    plt.close()
    # for T_idx, T in enumerate(THRESHOLDS):
    #     fig.add_trace(go.Scatter(x=MUTATIONS, y=sketch_recalls[T], mode="lines+markers+text", name=f"Sketch", textposition="top center", line=dict()))
    #     fig.add_trace(go.Scatter(x=MUTATIONS, y=graphaligner_recalls[T], mode="lines+markers+text", name=f"GraphAligner", textposition="top center", marker=dict(symbol=T_idx)))
    #     fig.add_trace(go.Scatter(x=MUTATIONS, y=vg_map_recalls[T], mode="lines+markers+text", name=f"vg map", textposition="top center", marker=dict(symbol=T_idx)))
    #     fig.add_trace(go.Scatter(x=MUTATIONS, y=vg_mpmap_recalls[T], mode="lines+markers+text", name=f"vg mpmap", textposition="top center", marker=dict(symbol=T_idx)))
    #     fig.update_layout(title="Baselines", xaxis_title="Mutation rate", yaxis_title=f"Recall @ T = {T}", font=dict(size=30))
    # fig.write_image(f"mutation_baselines_{T}.png", scale=1, width=1920, height=1080)
    print("Done")
