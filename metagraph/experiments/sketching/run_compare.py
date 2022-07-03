import plotly.graph_objects as go
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

seeder_map = {
    0: 'sketch',
    1: 'default',
    2: 'both'
}

BONUS = 5 # Change this in case we modify in metagraph
MATCH = 2
MISMATCH = 3
GAP_OPEN = 6
GAP_EXTENSION = 2

def get_metagraph_scores_and_time(result, num_seqs):
    scores = [0 for _ in range(num_seqs)]
    results = result.strip().split('\n')
    for i in range(num_seqs):
        res = json.loads(results[i])
        seq_num = int(res['name'][1:])
        if 'annotation' not in res:
            scores[seq_num] = -1e9
            continue
        c = list(Cigar(res['annotation']['cigar']).items())
        score = res['score']

        # Remove end bonuses
        if c[0][1] != 'S':
            score -= BONUS
        if c[-1][1] != 'S':
            score -= BONUS
        scores[seq_num] = score
    return np.asarray(scores), json.loads(results[-1])['time']


def get_gt_scores(references, mutations):
    assert len(references) == len(mutations)
    num = len(references)
    matrix = parasail.matrix_create("ACTG", MATCH, -MISMATCH)
    scores = []
    for i in range(0, num, 2):
        ref_seq = references[i+1].strip()
        mutated_seq = mutations[i+1].strip()
        score = parasail.nw(ref_seq, mutated_seq, GAP_OPEN, GAP_EXTENSION, matrix).score
        scores.append(score)
    return np.asarray(scores)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-k', type=int, required=True, default=15)
    parser.add_argument('--seeder', type=int, required=True)
    parser.add_argument('--embed-dim', type=int, default=20)
    parser.add_argument('--m', type=int, default=20)
    parser.add_argument('--stride', type=int, default=20)
    parser.add_argument('--n-times-sketch', type=int, default=1)
    parser.add_argument('--mutation-rate', type=int, default=1)
    parser.add_argument('--num-query-seqs', type=int, default=1000)
    parser.add_argument('--parallel', type=int, default=8)
    parser.add_argument('--batch-size', type=int, default=1000)
    parser.add_argument('--metagraph-path', type=str, required=True)
    parser.add_argument('--minimizer-window', type=int, required=True)
    parser.add_argument('--min-path-size', type=int, required=True)
    parser.add_argument('--max-path-size', type=int, required=True)
    parser.add_argument('--dataset-dir', type=str, required=True)

    args = parser.parse_args()

    if not os.path.exists("runs"):
        os.mkdir("runs")

    config = {
        'embed-dim': args.embed_dim,
        'n-times-sketch': args.n_times_sketch,
        'mutation-rate': args.mutation_rate,
        'num-query-seqs': args.num_query_seqs,
        'parallel': args.parallel,
        'batch-size': args.batch_size,
        'minimizer-window': args.minimizer_window,
        'min-path-size': args.min_path_size,
        'max-path-size': args.max_path_size
    }

    pprint(config, indent=4)


    K_VALS = range(20, args.max_k, 10)

    sketch_recalls = []
    sketch_times = []
    default_recalls = []
    default_times = []
    gt_scores = []
    if args.seeder == 0 or args.seeder == 2:
        print("Running for sketch...")
        # Sketch
        for K in tqdm(K_VALS):
            command = f"{args.metagraph_path} align " \
                          f"--seeder sketch " \
                          f"--output-path {args.dataset_dir} " \
                          f"--embed-dim {config['embed-dim']} " \
                          f"--n-times-sketch {config['n-times-sketch']} " \
                          f"--mutation-rate {config['mutation-rate']} " \
                          f"--minimizer-window {config['minimizer-window']} " \
                          f"--num-query-seqs {config['num-query-seqs']} " \
                          f"--min-path-size {config['min-path-size']} " \
                          f"--max-path-size {config['max-path-size']} " \
                          f"--parallel {config['parallel']} " \
                          f"--batch-size {config['batch-size']} " \
                          f"--json " \
                          f"--m {K // 2} " \
                          f"--stride {K // 4} " \
                          f"-i {args.dataset_dir}/sequence_{K}.dbg " \
                          "--experiment"
            print(command)
            result_ = subprocess.run(command.split(), capture_output=True, text=True)
            if len(gt_scores) == 0:
                gt_scores = get_gt_scores(open(args.dataset_dir + 'reference.fa', 'r').readlines(), open(args.dataset_dir + 'mutated.fa', 'r').readlines())

            metagraph_scores, time = get_metagraph_scores_and_time(result_.stdout, args.num_query_seqs)
            recall = np.sum(metagraph_scores >= gt_scores) / args.num_query_seqs
            sketch_recalls.append(recall)
            print(sketch_recalls, time)
            sketch_times.append(time)
    if args.seeder == 1 or args.seeder == 2:
        print("Running for default...")
        # Default
        for K in tqdm(K_VALS):
            command = f"{args.metagraph_path} align " \
                      f"--seeder default " \
                      f"--output-path {args.dataset_dir} " \
                      f"--embed-dim {config['embed-dim']} " \
                      f"--n-times-sketch {config['n-times-sketch']} " \
                      f"--mutation-rate {config['mutation-rate']} " \
                      f"--minimizer-window {config['minimizer-window']} " \
                      f"--num-query-seqs {config['num-query-seqs']} " \
                      f"--min-path-size {config['min-path-size']} " \
                      f"--max-path-size {config['max-path-size']} " \
                      f"--parallel {config['parallel']} " \
                      f"--batch-size {config['batch-size']} " \
                      f"--m {K // 2} " \
                      f"--stride {K // 4} " \
                      f"--json " \
                      f"-i {args.dataset_dir}/sequence_{K}.dbg " \
                      "--experiment"
            print(command)
            result_ = subprocess.run(command.split(), capture_output=True, text=True)
            if len(gt_scores) == 0:
                gt_scores = get_gt_scores(open(args.dataset_dir + 'reference.fa', 'r').readlines(), open(args.dataset_dir + 'mutated.fa', 'r').readlines())

            metagraph_scores, time = get_metagraph_scores_and_time(result_.stdout, args.num_query_seqs)
            print(time)
            recall = np.sum(metagraph_scores >= gt_scores) / args.num_query_seqs
            default_recalls.append(recall)
            default_times.append(time)
            print(default_recalls)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_times, y=sketch_recalls, text=list(K_VALS), textposition="top center", mode="lines+markers+text", name="Sketch"))
    fig.add_trace(go.Scatter(x=default_times, y=default_recalls, text=list(K_VALS), textposition="top center", mode="lines+markers+text", name="Default"))
    fig.update_xaxes(title_text="Average Time (s)", type="log")
    fig.update_yaxes(title_text="Recall")
    fig.add_annotation(text=json.dumps(config),
                       align='left',
                       showarrow=False,
                       xref='paper',
                       yref='paper',
                       x=0.0,
                       y=0.0,
                       bordercolor='black',
                       borderwidth=1)
    fig.write_image("recall.png", scale=1, width=1920, height=1080)