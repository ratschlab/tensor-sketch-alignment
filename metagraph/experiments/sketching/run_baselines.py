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

BONUS = 5 # Change this in case we modify in metagraph
MATCH = 2
MISMATCH = 3
GAP_OPEN = 6
GAP_EXTENSION = 2
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]#, 30]
# MUTATIONS = [0]#, 5, 10, 15, 20, 25]#, 30]
NUM_QUERY_SEQS = 300 


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
        cigar_obj = list(cigar.Cigar(base_cigar_string).items())
        if cigar_obj[0][1] == 'D':
            cigar_obj = cigar_obj[1:]
        print('base:', base_cigar_string)
        base_cigar_string = ''.join([str(x[0])+x[1] for x in cigar_obj])
        print('base post:', base_cigar_string)
        fixed_cigar_string = base_cigar_string.replace("S", "I")
        fixed_cigar_score = score_from_cigar(fixed_cigar_string)
        scores[seq_num] = fixed_cigar_score 
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

ALL_GT_SCORES = {}
if os.path.exists(DATASET_DIR + f'reference_0.fa'): #TODO: Make this nicer
    for mr in MUTATIONS:
        gt_scores, _ = get_gt_scores(open(DATASET_DIR+ f'reference_{mr}.fa', 'r').readlines(), open(DATASET_DIR + f'mutated_{mr}.fa', 'r').readlines())
        ALL_GT_SCORES[mr] = gt_scores
print(ALL_GT_SCORES)
def run_generate_seqs():
    print("Generating sequences...")
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} seqgen " \
                      f"--seeder sketch " \
                      f"--output-path {DATASET_DIR} " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {3*K + 10} " \
                      f"--max-path-size {3*K + 11} " \
                      f"--parallel 20 " \
                      f"--batch-size 1000 " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      "--experiment"
        print(command)
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)

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
                      f"--mutation-rate {mr} " \
                      f"--parallel 20 " \
                      f"--batch-size 1000 " \
                      f"--json " \
                      f"--m {K} " \
                      f"--stride 20 " \
                      f"--num-neighbours 10 " \
                      f"--align-end-bonus {BONUS} " \
                      f"--align-match-score {MATCH} " \
                      f"--align-mm-transition-penalty {MISMATCH} " \
                      f"--align-mm-transversion-penalty {MISMATCH} " \
                      f"--align-gap-open-penalty {GAP_OPEN} " \
                      f"--align-gap-extension-penalty {GAP_EXTENSION} " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        print(command)
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        end = mytime.time()
        gt_scores = ALL_GT_SCORES[mr]
        metagraph_scores, time, precision = get_metagraph_scores_and_time_and_precision(result_.stdout, NUM_QUERY_SEQS)
        print(gt_scores, metagraph_scores)
        recall = np.sum(metagraph_scores >= gt_scores) / NUM_QUERY_SEQS 
        sketch_recalls.append(recall)
        print(sketch_recalls)
        sketch_times.append(time/NUM_QUERY_SEQS)
    return sketch_recalls, sketch_times

def run_metagraph_default():
    print("Running default...")
    sketch_recalls = []
    sketch_times = []
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} align " \
                      f"--seeder default " \
                      f"--output-path {DATASET_DIR} " \
                      f"--embed-dim 14 " \
                      f"--n-times-sketch 1 " \
                      f"--mutation-rate {mr} " \
                      f"--parallel 20 " \
                      f"--batch-size 1000 " \
                      f"--json " \
                      f"--m {K} " \
                      f"--stride 20 " \
                      f"--num-neighbours 10 " \
                      f"--align-end-bonus {BONUS} " \
                      f"--align-match-score {MATCH} " \
                      f"--align-mm-transition-penalty {MISMATCH} " \
                      f"--align-mm-transversion-penalty {MISMATCH} " \
                      f"--align-gap-open-penalty {GAP_OPEN} " \
                      f"--align-gap-extension-penalty {GAP_EXTENSION} " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      f"{DATASET_DIR}/mutated_{mr}.fa "
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)
        end = mytime.time()
        gt_scores = ALL_GT_SCORES[mr]
        metagraph_scores, time, precision = get_metagraph_scores_and_time_and_precision(result_.stdout, NUM_QUERY_SEQS)
        print(metagraph_scores, gt_scores)
        recall = np.sum(metagraph_scores >= gt_scores) / NUM_QUERY_SEQS 
        sketch_recalls.append(recall)
        sketch_times.append(time/NUM_QUERY_SEQS)
    return sketch_recalls, sketch_times

import re
def convert_to_cigar(s):
    # - deletion
    # + insertion
    # :3 matches
    # * substitution
    cigar_string = ""
    regex_expr = "-[AGCT]+|\+[AGCT]+|:[0-9]+|\*[ACGT]+"
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

def run_vg_map():
    print("Running vg mapping...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        f = open(f'data/out_{mr}.gam','w')
        scores = np.zeros(NUM_QUERY_SEQS) - np.inf
        gt_scores = ALL_GT_SCORES[mr]
        # command = f"{VG_PATH} map --full-l-bonus 0 --match {MATCH} --mismatch {MISMATCH} --gap-open {GAP_OPEN} --gap-extend {GAP_EXTENSION} -f data/mutated_{mr}.fa -x data/sequence.xg -g sequence.gcsa --gaf"
        # command = f"{VG_PATH} map -f data/mutated_{mr}.fa -x data/sequence.xg -g sequence.gcsa --gaf"
        # command = f"{VG_PATH} map --drop-full-l-bonus --match {MATCH} --mismatch {MISMATCH} --gap-open {GAP_OPEN} --gap-extend {GAP_EXTENSION} -f data/mutated_{mr}.fa -x data/sequence.xg -g data/sequence.gcsa"
        # ./vg mpmap -n DNA -F GAF -x data/sequence.xg -g data/sequence.gcsa -f data/mutated_20.fa
        command = f"{VG_PATH} mpmap --match {MATCH} --full-l-bonus {BONUS} --mismatch {MISMATCH} --gap-open {GAP_OPEN} --gap-extend {GAP_EXTENSION} -n DNA -F GAF -x data/sequence.xg -g data/sequence.gcsa -f data/mutated_{mr}.fa"
        start = mytime.time()
        subprocess.run(command.split(), stdout=f)
        end = mytime.time()
        f.close()
        # command = f"{VG_PATH} view -a data/out_{mr}.gam"
        # alignments = subprocess.run(command.split(), capture_output=True).stdout.decode('utf-8').replace('\n',',')
        # alignments = json.loads('[' + alignments[:-1] + ']')
        # print(alignments)
        # scores = parse_scores(alignments)
        alignments = open(f"data/out_{mr}.gam", "r").readlines()
        for alignment in alignments:
            parts = alignment.split("\t")
            if len(parts) < 13:
                continue
            query = int(parts[0][1:])
            csz_string = parts[13][5:]
            cigar_str = convert_to_cigar(csz_string)
            score = score_from_cigar(cigar_str)
            scores[query] = score
        scores = np.asarray(scores)
        recall = np.sum(scores >= gt_scores) / NUM_QUERY_SEQS 
        print(gt_scores, scores)
        recalls.append(recall)
        print(recalls)
        times.append((end-start)/NUM_QUERY_SEQS)
    return np.array(recalls), np.array(times)

def run_graphaligner_vg():
    print("Running graphaligner vg...")
    recalls = []
    times = []
    for mr in tqdm(MUTATIONS):
        gt_scores = ALL_GT_SCORES[mr]
        command = f"GraphAligner --try-all-seeds --seeds-mem-count -1 -g data/sequence_{K}_blunted.gfa -C -1 -f data/mutated_{mr}.fa -a data/aln_{mr}.gaf -b 10000" #-x dbg"
        start = mytime.time()
        subprocess.run(command.split())
        end = mytime.time()
        
        scores = np.zeros(NUM_QUERY_SEQS) - np.inf

        alignments = csv.reader(open(f"data/aln_{mr}.gaf", 'r'), delimiter='\t')
        for alignment in alignments:
            cigar_str = alignment[-1][5:]
            score = score_from_cigar(cigar_str)
            idx = int(alignment[0][1:])
            scores[idx] = score

        recall = np.sum(scores >= gt_scores) / NUM_QUERY_SEQS 
        recalls.append(recall)
        times.append(end-start)
    return np.array(recalls), np.array(times)

def run_minigraph():
    print("Running minigraph...")
    recalls = []
    times = []

    command = f"{MINIGRAPH_PATH} -cxggs -t16 data/sequence.fa"
    f = open("data/sequence.gfa", 'w')
    subprocess.run(command.split(), stdout=f)
    f.close()

    #minigraph -cxggs -t16 ref.fa sample1.fa sample2.fa > out.gfa
    for mr in tqdm(MUTATIONS):
        gt_scores = ALL_GT_SCORES[mr]
        # gt_scores, _ = get_gt_scores(open(dataset_dir + f'reference_{mr}.fa', 'r').readlines(), open(dataset_dir + f'mutated_{mr}.fa', 'r').readlines())
        command = f"{MINIGRAPH_PATH} -c -k 8 -w 1 data/sequence.gfa data/mutated_{mr}.fa"
        start = mytime.time()
        subprocess.run(command.split(), stdout = open(f"data/minigraph_{mr}.gaf", 'w'))
        end = mytime.time()

        scores = np.zeros(NUM_QUERY_SEQS) - np.inf
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
    # default_recalls, default_times = run_metagraph_default()
    vg_recalls, vg_times = run_vg_map()
    # ga_recalls, ga_times = run_graphaligner_vg()
    # mini_recalls, mini_times = run_minigraph()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_times, y=sketch_recalls, mode="lines+markers+text", name="Sketch", text=MUTATIONS, textposition="top center"))
    # fig.add_trace(go.Scatter(x=default_times, y=default_recalls, mode="lines+markers+text", name="Default metagraph", text=MUTATIONS, textposition="top center"))
    fig.add_trace(go.Scatter(x=vg_times, y=vg_recalls, mode="lines+markers+text", name="vg mpmap", text=MUTATIONS, textposition="top center"))
    # # fig.add_trace(go.Scatter(x=ga_times, y=ga_recalls, mode="lines+markers+text", name="GraphAligner", text=MUTATIONS, textposition="top center"))
    # # fig.add_trace(go.Scatter(x=mini_times, y=mini_recalls, mode="lines+markers+text", name="minigraph", text=MUTATIONS, textposition="top center"))
    fig.update_layout(title="Baselines", xaxis_title="Average time(s)", yaxis_title="Recall @ mutation rates", font=dict(size=25))
    fig.write_image("baselines.png", scale=1, width=1920, height=1080)
    print("Done")
