import cigar
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


BONUS = 5 # Change this in case we modify in metagraph
MATCH = 2
MISMATCH = 3
GAP_OPEN = 6
GAP_EXTENSION = 2
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25, 30]
NUM_QUERY_SEQS = 100

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

def get_metagraph_scores_and_time_and_precision(result, num_seqs):
    scores = [0 for _ in range(num_seqs)]
    results = result.strip().split('\n')
    for i in range(num_seqs):
        res = json.loads(results[i])
        seq_num = int(res['name'][1:])
        if 'annotation' not in res:
            scores[seq_num] = -1
            continue
        c = list(Cigar(res['annotation']['cigar']).items())
        score = res['score']
        # Remove end bonuses
        if c[0][1] != 'S':
            score -= BONUS
        if c[-1][1] != 'S':
            score -= BONUS
        scores[seq_num] = score
    precision = json.loads(results[-1])['precision']
    return np.asarray(scores), json.loads(results[-1])['time'], precision

def get_gt_scores(references, mutations):
    assert len(references) == len(mutations)
    num = len(references)
    matrix = parasail.matrix_create("ACTG", MATCH, -MISMATCH)
    scores = []
    dists = []
    for i in range(0, num, 2):
        ref_seq = references[i+1].strip()
        mutated_seq = mutations[i+1].strip()
        score = parasail.nw(ref_seq, mutated_seq, GAP_OPEN, GAP_EXTENSION, matrix).score
        dists.append(edlib.align(ref_seq, mutated_seq, task='path')['editDistance'])
        scores.append(score)
    return np.asarray(scores), np.asarray(dists)

def run_generate_seqs():
    print("Generating sequences...")
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder sketch " \
                      f"--output-path {dataset_dir} " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {2*K + 10} " \
                      f"--max-path-size {2*K + 11} " \
                      f"--parallel 20 " \
                      f"--batch-size 20 " \
                      f"--json " \
                      f"--m {K} " \
                      f"--stride 20 " \
                      f"--num-neighbours 10 " \
                      f"-i {dataset_dir}/sequence_{K}.dbg " \
                      "--experiment"
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)

def run_metagraph_sketching():
    print("Running Sketching...")
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    sketch_recalls = []
    sketch_times = []
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder sketch " \
                      f"--output-path {dataset_dir} " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {2*K + 10} " \
                      f"--max-path-size {2*K + 11} " \
                      f"--parallel 20 " \
                      f"--batch-size 20 " \
                      f"--json " \
                      f"--m {K} " \
                      f"--stride 20 " \
                      f"--num-neighbours 10 " \
                      f"-i {dataset_dir}/sequence_{K}.dbg " \
                      f"{dataset_dir}/mutated_{mr}.fa "
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        end = mytime.time()
        gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        metagraph_scores, time, precision = get_metagraph_scores_and_time_and_precision(result_.stdout, NUM_QUERY_SEQS)
        recall = np.sum(metagraph_scores >= gt_scores) / NUM_QUERY_SEQS 
        sketch_recalls.append(recall)
        sketch_times.append(time)
    return sketch_recalls, sketch_times

def run_metagraph_default():
    print("Running default...")
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    sketch_recalls = []
    sketch_times = []
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder default " \
                      f"--output-path {dataset_dir} " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {2*K + 10} " \
                      f"--max-path-size {2*K + 11} " \
                      f"--parallel 20 " \
                      f"--batch-size 20 " \
                      f"--json " \
                      f"--m {K} " \
                      f"--stride 20 " \
                      f"--num-neighbours 10 " \
                      f"-i {dataset_dir}/sequence_{K}.dbg " \
                      f"{dataset_dir}/mutated_{mr}.fa "
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        end = mytime.time()
        gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        metagraph_scores, time, precision = get_metagraph_scores_and_time_and_precision(result_.stdout, NUM_QUERY_SEQS)
        recall = np.sum(metagraph_scores >= gt_scores) / NUM_QUERY_SEQS 
        sketch_recalls.append(recall)
        sketch_times.append(time)
    return sketch_recalls, sketch_times

def parse_scores(alignments):
    scores = np.zeros(NUM_QUERY_SEQS)
    for alignment in alignments:
        idx = int(alignment['name'][1:])

        if 'score' in alignment:
            score = alignment['score']
            scores[idx] = score
    return scores

def run_vg_map():
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    print("Running vg mapping...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        f = open(f'data/out_{mr}.gam','w')
        gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        command = f"{VG_PATH} map --full-l-bonus 0 --match {MATCH} --mismatch {MISMATCH} --gap-open {GAP_OPEN} --gap-extend {GAP_EXTENSION} -f data/mutated_{mr}.fa -x data/sequence.xg -g sequence.gcsa"
        start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        end = mytime.time()
        f.close()
        command = f"{VG_PATH} view -a data/out_{mr}.gam"
        alignments = subprocess.run(command.split(), capture_output=True).stdout.decode('utf-8').replace('\n',',')
        alignments = json.loads('[' + alignments[:-1] + ']')
        scores = parse_scores(alignments)
        recall = np.sum(scores >= gt_scores) / NUM_QUERY_SEQS 
        print(scores, gt_scores)
        recalls.append(recall)
        times.append(end-start)
    return np.array(recalls), np.array(times)

# GraphAligner -g sequence_50.gfa -f mutated_10.fa -a aln.gaf -x dbg
import csv
def run_graphaligner_vg():
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    print("Running graphaligner vg...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        command = f"GraphAligner --try-all-seeds -g data/sequence_{K}.gfa -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf --seeds-mxm-length 9 -b 30 --seeds-mum-count -1" #-x dbg"
        start = mytime.time()
        subprocess.run(command.split())
        end = mytime.time()
        
        scores = np.zeros(NUM_QUERY_SEQS)

        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            cigar_str = alignment[-1][5:]
            score = score_from_cigar(cigar_str)
            idx = int(alignment[0][1:])
            scores[idx] = score

        print(scores, gt_scores)
        recall = np.sum(scores >= gt_scores) / NUM_QUERY_SEQS 
        recalls.append(recall)
        times.append(end-start)
    return np.array(recalls), np.array(times)

def run_minigraph():
    dataset_dir = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
    print("Running minigraph...")
    recalls = []
    times = []

    command = f"{MINIGRAPH_PATH} -cxggs -t16 data/sequence.fa"
    f = open("data/sequence.gfa", 'w')
    subprocess.run(command.split(), stdout=f)
    f.close()

    #minigraph -cxggs -t16 ref.fa sample1.fa sample2.fa > out.gfa
    for mr in tqdm(MUTATIONS):
        gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        command = f"{MINIGRAPH_PATH} -c -k 8 -w 1 data/sequence.gfa data/mutated_{mr}.fa"
        start = mytime.time()
        subprocess.run(command.split(), stdout = open(f"data/minigraph_{mr}.gaf", 'w'))
        end = mytime.time()

        scores = np.zeros(NUM_QUERY_SEQS)
        alignments = csv.reader(open(f"data/minigraph_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            cigar_str = alignment[-1][5:]
            score = score_from_cigar(cigar_str)
            idx = int(alignment[0][1:])
            scores[idx] = score
        recall = np.sum(scores >= gt_scores) / NUM_QUERY_SEQS 
        recalls.append(recall)
        times.append(end-start)
    return np.array(recalls), np.array(times)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", action='store_true')
    args = parser.parse_args()
    if args.gen:
        run_generate_seqs()
        exit(0)
    sketch_recalls, sketch_times = run_metagraph_sketching()
    default_recalls, default_times = run_metagraph_default()
    vg_recalls, vg_times = run_vg_map()
    ga_recalls, ga_times = run_graphaligner_vg()
    mini_recalls, mini_times = run_minigraph()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_times, y=sketch_recalls, mode="lines+markers+text", name="Sketch", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=default_times, y=default_recalls, mode="lines+markers+text", name="Default metagraph", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=vg_times, y=vg_recalls, mode="lines+markers+text", name="vg map", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=ga_times, y=ga_recalls, mode="lines+markers+text", name="GraphAligner", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=mini_times, y=mini_recalls, mode="lines+markers+text", name="minigraph", text=MUTATIONS, textposition="top center"))
    fig.update_layout(title="Baselines", xaxis_title="Cumulative time(s)", yaxis_title="Recall @ mutation rates")
    fig.write_image("baselines.png", scale=1, width=1920, height=1080)
    print("Done")
