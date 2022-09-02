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
import edlib

METAGRAPH_PATH = "/home/alex/metagraph/metagraph/build/metagraph"
ASTARIX_PATH = "/home/alex/benchmark/datagen/astarix/release/astarix"
VG_PATH = "/home/alex/benchmark/datagen/vg"
MINIGRAPH_PATH = "/home/alex/benchmark/datagen/minigraph/minigraph"
DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"

BONUS = 0 # Change this in case we modify in metagraph
MATCH = 0
MISMATCH = 1
GAP_OPEN = 1
GAP_EXTENSION = 1
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]
# MUTATIONS = [5]#, 10, 15, 20, 25]
NUM_QUERY_SEQS = 300


def get_gt_scores(references, mutations):
    assert len(references) == len(mutations)
    num = len(references)
    matrix = parasail.matrix_create("ACTG", MATCH, -MISMATCH)
    scores = []
    dists = []
    for i in range(0, num, 2):
        ref_seq = references[i+1].strip()
        mutated_seq = mutations[i+1].strip()
        alignment = parasail.nw_trace(ref_seq, mutated_seq, GAP_OPEN, GAP_EXTENSION, matrix)
        score = alignment.score
        dists.append(edlib.align(ref_seq, mutated_seq, task='path')['editDistance'])
        scores.append(score)
    return np.asarray(scores), np.asarray(dists)

ALL_GT_SCORES= {}
for mr in MUTATIONS:
    gt_scores, gt_dists = get_gt_scores(open(DATASET_DIR+ f'reference_{mr}.fa', 'r').readlines(), open(DATASET_DIR + f'mutated_{mr}.fa', 'r').readlines())
    ALL_GT_SCORES[mr] = gt_scores

def score_from_cigar(s):
    score = 0
    parts = list(cigar.Cigar(s).items())
    n_parts = len(parts)
    i = 0
    while i < n_parts:
        count, op = parts[i]
        if op == '=':
            score += count * MATCH
        elif op == 'X':
            score -= count * MISMATCH
        elif op == 'I' or op == 'D':
            score -= (GAP_OPEN + (count-1) * GAP_EXTENSION)
        i += 1
    return score

def parse_scores(alignments):
    scores = np.zeros(NUM_QUERY_SEQS) - np.inf
    for alignment in alignments:
        idx = int(alignment['name'][1:])
        if 'score' in alignment:
            score = alignment['score']
            scores[idx] = score
    return scores

def get_metagraph_scores_and_time_and_precision(result, num_seqs):
    scores = np.zeros(num_seqs) - np.inf 
    results = result.strip().split('\n')
    for i in range(num_seqs):
        res = json.loads(results[i])
        seq_num = int(res['name'][1:])
        if 'annotation' not in res:
            continue
        base_cigar_string = res['annotation']['cigar']
        print("Q" + str(seq_num))
        print("base: ", base_cigar_string)
        fixed_cigar_string = base_cigar_string.replace("S", "I")
        print("fixed: ", fixed_cigar_string)
        fixed_cigar_score = score_from_cigar(fixed_cigar_string)
        scores[seq_num] = fixed_cigar_score 
    precision = json.loads(results[-1])['precision']
    return np.asarray(scores), json.loads(results[-1])['time'], precision


def run_metagraph_sketching():
    print("Running Sketching...")
    sketch_recalls = []
    sketch_times = []
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder sketch " \
                      f"--output-path {DATASET_DIR} " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--parallel 20 " \
                      f"--batch-size 1000 " \
                      f"--num-neighbours 10 " \
                      f"--json " \
                      f"--align-edit-distance " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        print(command)
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        print(result)
        end = mytime.time()
        gt_score = ALL_GT_SCORES[mr]
        metagraph_scores, time, precision = get_metagraph_scores_and_time_and_precision(result_.stdout, NUM_QUERY_SEQS)
        print(metagraph_scores, gt_score)
        recall = np.sum(metagraph_scores >= gt_score) / NUM_QUERY_SEQS 
        sketch_recalls.append(recall)
        print(sketch_recalls)
        sketch_times.append(time/NUM_QUERY_SEQS)
    return sketch_recalls, sketch_times

def run_metagraph_sketching2():
    print("Running Sketching...")
    sketch_recalls = []
    sketch_times = []
    for mr in tqdm(MUTATIONS):
        ref_seqs = open(DATASET_DIR+ f'reference_{mr}.fa', 'r').readlines()
        query_seqs = open(DATASET_DIR+ f'mutated_{mr}.fa', 'r').readlines()
        scores = init_scores(query_seqs) 
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder sketch " \
                      f"--output-path {DATASET_DIR} " \
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
            # dists.append(edlib.align(ref_seq, mutated_seq, task='path')['editDistance'])
            seq_idx = int(alignment[0][1:])
            matched_seq = alignment[3]
            ref_seq = ref_seqs[seq_idx * 2 + 1].strip()
            scores[seq_idx] = edlib.align(ref_seq, matched_seq, task='path')['editDistance']
        recall = np.mean(scores)
        sketch_recalls.append(recall)
        sketch_times.append(time/NUM_QUERY_SEQS)
    return sketch_recalls, sketch_times


def run_graphaligner():
    print("Running graphaligner...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        gt_score = ALL_GT_SCORES[mr]
        command = f"GraphAligner --try-all-seeds --seeds-mem-count -1 -g data/sequence_{K}_blunted.gfa -C -1 -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -b 10000" #-x dbg"
        start = mytime.time()
        subprocess.run(command.split())
        end = mytime.time()
        scores = np.zeros(NUM_QUERY_SEQS) - np.inf
        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        print("GRAPH")
        for alignment in alignments:
            print(alignment)
            break
            cigar_str = alignment[-1][5:]
            print(cigar_str)
            score = score_from_cigar(cigar_str)
            idx = int(alignment[0][1:])
            scores[idx] = score
        recall = np.sum(scores >= gt_score) / NUM_QUERY_SEQS 
        recalls.append(recall)
        times.append((end-start) / NUM_QUERY_SEQS)
    return np.array(recalls), np.array(times)

def init_scores(query_seqs):
    scores = []
    for i in range(0, len(query_seqs), 2):
        header = query_seqs[i].strip()
        seq = query_seqs[i+1].strip()
        scores.append(len(seq))
    return np.array(scores)

def run_graphaligner2():
    print("Running graphaligner...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        query_seqs = open(DATASET_DIR+ f'mutated_{mr}.fa', 'r').readlines()
        scores = init_scores(query_seqs) 
        gt_score = ALL_GT_SCORES[mr]
        command = f"GraphAligner --try-all-seeds --seeds-mem-count -1 -g data/sequence_{K}_blunted.gfa -C -1 -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -b 10000" #-x dbg"
        start = mytime.time()
        subprocess.run(command.split())
        end = mytime.time()
        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            seq_idx = int(alignment[0][1:])
            nm = float(alignment[12][5:])
            scores[seq_idx] = nm
        recall = np.mean(scores) 
        recalls.append(recall)
        print(recalls)
        times.append((end-start) / NUM_QUERY_SEQS)
    return np.array(recalls), np.array(times)

if __name__ == '__main__':
    sketch_recalls, sketch_times = run_metagraph_sketching2()
    graphaligner_recalls, graphaligner_times = run_graphaligner2()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_times, y=sketch_recalls, mode="lines+markers+text", name="Sketch", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=graphaligner_times, y=graphaligner_recalls, mode="lines+markers+text", name="GraphAligner", text=MUTATIONS, textposition="top center"))
    fig.update_layout(title="Baselines", xaxis_title="Average time(s)", yaxis_title="Recall @ mutation rates", font=dict(size=30))
    fig.write_image("edit_distance_baselines.png", scale=1, width=1920, height=1080)
    print("Done")
