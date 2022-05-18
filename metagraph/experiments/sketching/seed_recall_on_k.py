import os
import numpy as np
import subprocess
import random
import plotly
import json
import plotly.graph_objects as go

DATASET_DIR = './data'
METAGRAPH_PATH = "/Users/alex/metagraph/metagraph/build/metagraph"
MAX_K = 80

if __name__ == '__main__':
    x = []
    y = []
    for K in range(10, MAX_K, 10):
        print(K)
        command = f"{METAGRAPH_PATH} align " \
                  "--seeder sketch " \
                  "--sketch_dim 10 " \
                  "--n_times_subsample 10 " \
                  "--subsampled_sketch_dim 8 " \
                  "--mutation-rate 20 " \
                  "-v " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  f"{DATASET_DIR}/query_strings.fa"
        result = subprocess.run(command.split(), capture_output=True, text=True)
        output = json.loads(result.stdout.strip().split('\n')[-1])
        x.append(output['avg_time'])
        y.append(output['recall'])
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y, name=f"Recall")
    )
    fig.update_layout(
        title="Seed recall",
        xaxis_title="Average time",
        yaxis_title="Recall",
    )
    fig.show()
