import os
import numpy as np
import subprocess
import random
import plotly
import json
import plotly.graph_objects as go

DATASET_DIR = './data'
METAGRAPH_PATH = "/Users/alex/metagraph/metagraph/build/metagraph"
MAX_K = 84


def recall(result):
    n_alignments = len(result)
    n_correct_aligned = 0
    for alignment in result:
        #print(alignment)
        alignment = json.loads(alignment)
        if "read_mapped" in alignment and alignment['read_mapped'] is True:
            n_correct_aligned += 1
    return n_correct_aligned / n_alignments


if __name__ == '__main__':
    x = []
    y = []
    for K in range(2, MAX_K, 1):
        print(K)
        command = f"{METAGRAPH_PATH} align " \
                  "--seeder sketch " \
                  "--sketch_dim 5 " \
                  "--n_times_subsample 5 " \
                  "--subsampled_sketch_dim 3 " \
                  "--json " \
                  "-v " \
                  f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                  f"{DATASET_DIR}/query_strings.fa"
        result = subprocess.run(command.split(), capture_output=True, text=True)
        print(result.stdout)
        recall_result = recall(result.stdout.strip().split("\n"))
        x.append(K)
        y.append(recall_result)
        print(recall_result)
        #exit(0)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y, name=f"Recall")
    )
    fig.show()
