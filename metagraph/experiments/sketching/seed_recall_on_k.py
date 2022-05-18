import subprocess
import json
import plotly.graph_objects as go
import time
import os
from time import localtime, strftime

DATASET_DIR = './data'
METAGRAPH_PATH = "/Users/alex/metagraph/metagraph/build/metagraph"
MAX_K = 90

if __name__ == '__main__':
    if not os.path.exists("runs"):
        os.mkdir("runs")

    x = []
    y = []
    K_VALS = list(range(11, MAX_K, 10))
    total_time = 0

    config = {
        'sketch_dim': 15,
        'n_times_subsample': 2,
        'subsampled_sketch_dim': 13,
        'mutation-rate': 40,
        'num-query-seqs': 1000
    }

    for K in K_VALS:
        start = time.time()
        print(K)
        command = f"{METAGRAPH_PATH} align " \
                  "--seeder sketch " \
                  f"--sketch_dim {config['sketch_dim']} " \
                  f"--n_times_subsample {config['n_times_subsample']} " \
                  f"--subsampled_sketch_dim {config['subsampled_sketch_dim']} " \
                  f"--mutation-rate {config['mutation-rate']} " \
                  f"--num-query-seqs {config['num-query-seqs']} " \
                  "-v " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  "--experiment"
        result = subprocess.run(command.split(), capture_output=True, text=True)
        output = json.loads(result.stdout.strip().split('\n')[-1])
        x.append(output['avg_time'])
        y.append(output['recall'])
        end = time.time()

        print(f"Time: {(end - start):.2f}\n")

        total_time += (end - start)

    print(f"Experiment total time: {total_time:.2f}")

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y, text=K_VALS, textposition="top center", mode="lines+markers+text", name=f"Recall")
    )
    fig.update_layout(
        title="Seed recall",
        xaxis_title="Average time",
        yaxis_title="Recall",
    )
    fig.show()

    experiment_dir = os.path.join("runs", strftime("%Y-%m-%dT%H:%M:%S", localtime()))
    os.mkdir(experiment_dir)

    fig.write_image(os.path.join(experiment_dir, 'fig.png'), scale=1, width=1920, height=1080)
    with open(os.path.join(experiment_dir, 'config.json'), 'w') as outfile:
        json.dump(config, outfile, indent=4)

    print(f"Logged experiment at: {experiment_dir}")
