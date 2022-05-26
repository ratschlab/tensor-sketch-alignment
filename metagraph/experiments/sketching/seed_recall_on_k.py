import subprocess
import json
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import time
import os
from time import localtime, strftime
import argparse
from pprint import pprint

DATASET_DIR = './data'
LOGS = []


def print2(s):
    print(s)
    LOGS.append(str(s))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--embed-dim', type=int, default=30)
    parser.add_argument('--n-times-sketch', type=int, default=20)
    parser.add_argument('--mutation_rate', type=int, default=15)
    parser.add_argument('--num_query_seqs', type=int, default=1000)
    parser.add_argument('--parallel', type=int, default=8)
    parser.add_argument('--batch-size', type=int, default=500)
    parser.add_argument('--seeder', type=str, default="sketch")
    parser.add_argument('--metagraph-path', type=str, required=True)
    parser.add_argument('--max-k', type=int, required=True)
    parser.add_argument('--output-path', type=str, required=True)
    parser.add_argument('--min-path-size', type=int, required=True)
    parser.add_argument('--max-path-size', type=int, required=True)
    args = parser.parse_args()

    METAGRAPH_PATH = args.metagraph_path
    MAX_K = args.max_k

    if not os.path.exists("runs"):
        os.mkdir("runs")

    x = []
    recall = []
    precision = []
    K_VALS = list(range(11, MAX_K, 10))
    total_time = 0

    config = {
        'embed-dim': args.embed_dim,
        'n-times-sketch': args.n_times_sketch,
        'mutation-rate': args.mutation_rate,
        'num-query-seqs': args.num_query_seqs,
        'parallel': args.parallel,
        'batch-size': args.batch_size,
        'seeder': args.seeder,
        'min-path-size': args.min_path_size,
        'max-path-size': args.max_path_size
    }

    print2("Launching experiment")
    pprint(config)

    for K in K_VALS:
        start = time.time()
        print2(K)
        command = f"{METAGRAPH_PATH} align " \
                  f"--seeder {config['seeder']} " \
                  f"--output-path {args.output_path} " \
                  f"--embed-dim {config['embed-dim']} " \
                  f"--n-times-sketch {config['n-times-sketch']} " \
                  f"--mutation-rate {config['mutation-rate']} " \
                  f"--num-query-seqs {config['num-query-seqs']} " \
                  f"--min-path-size {config['min-path-size']} " \
                  f"--max-path-size {config['max-path-size']} " \
                  f"--parallel {config['parallel']} " \
                  f"--batch-size {config['batch-size']} " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  "--experiment"
        print2(command)
        result = subprocess.run(command.split(), capture_output=True, text=True)
        output = json.loads(result.stdout.strip().split('\n')[-1])
        x.append(output['avg_time'])
        recall.append(output['recall'])
        precision.append(output['precision'])
        print2(x)
        print2(recall)
        print2(precision)
        end = time.time()

        print2(f"Time: {(end - start):.2f}\n")

        total_time += (end - start)

    print2(f"Experiment total time: {total_time:.2f}")

    fig = make_subplots(rows=2, cols=1)
    fig.add_trace(
        go.Scatter(x=x, y=recall, text=K_VALS, textposition="top center", mode="lines+markers+text", name=f"Recall"),
        row=1,
        col=1
    )
    fig.add_trace(
        go.Scatter(x=x, y=precision, text=K_VALS, textposition="top center", mode="lines+markers+text", name=f"Precision"),
        row=2,
        col=1
    )

    fig.update_xaxes(title_text="Average Time (s)", row=1, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=2, col=1)
    fig.update_yaxes(title_text="Recall", row=1, col=1)
    fig.update_yaxes(title_text="Precision", row=3, col=1)

    experiment_dir = os.path.join("runs", strftime("%Y-%m-%dT%H:%M:%S", localtime()))
    os.mkdir(experiment_dir)

    fig.write_image(os.path.join(experiment_dir, 'fig.png'), scale=1, width=1920, height=1080)
    with open(os.path.join(experiment_dir, 'config.json'), 'w') as outfile:
        json.dump(config, outfile, indent=4)
    data = {
        'avg_time': x,
        'recall': recall,
        'precision': precision,
    }

    with open(os.path.join(experiment_dir, 'points.json'), 'w') as outfile:
        json.dump(data, outfile, indent=4)

    print2(f"Logged experiment at: {experiment_dir}")
    with open(os.path.join(experiment_dir, 'logs.txt'), 'w') as f:
        f.write('\n'.join(LOGS))
