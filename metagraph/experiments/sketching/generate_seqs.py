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
DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]
NUM_QUERY_SEQS = 500 

def run_generate_seqs():
    print("Generating sequences...")
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} seqgen " \
                      f"--output-path {DATASET_DIR} " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {3*K + 10} " \
                      f"--max-path-size {3*K + 11} " \
                      f"-i {DATASET_DIR}/sequence_{K}.dbg " \
                      "--experiment"
        print(command)
        start = mytime.time()
        result_ = subprocess.run(command.split(), capture_output=True, text=True)

if __name__ == '__main__':
    run_generate_seqs()
    print(f"Num seqs: {NUM_QUERY_SEQS}")
    print("Done")
