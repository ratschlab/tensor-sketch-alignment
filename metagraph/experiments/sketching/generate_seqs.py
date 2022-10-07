import argparse
import os
import subprocess
import numpy as np
from tqdm import tqdm

#METAGRAPH_PATH = "/home/alex/metagraph/metagraph/build/metagraph"
#DATASET_DIR = "/home/alex/metagraph/metagraph/experiments/sketching/data/"

METAGRAPH_PATH = "/cluster/apps/biomed/grlab/ameterez/metagraph/metagraph/build/metagraph"
DATASET_DIR = None
K = 80
MUTATIONS = [0, 5, 10, 15, 20, 25]
NUM_QUERY_SEQS = 500 

def run_generate_seqs():
    save_path = os.path.join(DATASET_DIR, f"sequence_{K}.dbg")
    
    print("Generating sequences...")
    for mr in tqdm(MUTATIONS):
        command = f"{METAGRAPH_PATH} seqgen " \
                      f"--output-path {DATASET_DIR} " \
                      f"--mutation-rate {mr} " \
                      f"--num-query-seqs {NUM_QUERY_SEQS} " \
                      f"--min-path-size {3*K + 10} " \
                      f"--max-path-size {3*K + 11} " \
                      f"-i {save_path} " \
                      "--experiment"
        print(command)
        result_ = subprocess.run(command.split(), capture_output=True, text=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", type=str, required=True)
    args = parser.parse_args()

    DATASET_DIR = args.dataset_dir + "/"
    run_generate_seqs()
    print(f"Num seqs: {NUM_QUERY_SEQS}")
    print("Done")
