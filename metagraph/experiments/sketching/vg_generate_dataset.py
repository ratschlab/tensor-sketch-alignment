import os
import numpy as np
import subprocess
import random
import shutil
import argparse


DATASET_DIR = None 
get_blunted_path = "/cluster/apps/biomed/grlab/ameterez/GetBlunted/build/get_blunted"
vg_path = "/cluster/apps/biomed/grlab/ameterez/vg"

INPUT_SEQ = "sequence.fa"
K = 80

if __name__ == '__main__':
    assert os.path.exists(vg_path), print(vg_path)
    assert os.path.exists(get_blunted_path), print(get_blunted_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", type=str, required=True)
    args = parser.parse_args()
    DATASET_DIR = args.dataset_dir
    print(DATASET_DIR)
    assert os.path.exists(DATASET_DIR), "Please create dataset directory"
    
    dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}')
    blunted_dbg_output = os.path.join(DATASET_DIR, INPUT_SEQ.split('.')[0] + f'_{K}_blunted')

    # Blunt graph
    blunt_command = f"{get_blunted_path} -t 8 --input_gfa {dbg_output}.gfa"
    blunted_graph = subprocess.run(blunt_command.split(), capture_output=True).stdout.decode("utf-8")
    open(f"{blunted_dbg_output}.gfa", 'w').write(blunted_graph)
   
    # Build xg index
    vg_command = f"{vg_path} autoindex -t 8 -T /cluster/work/grlab/ameterez/temp/ -g {DATASET_DIR}/sequence_{K}_blunted.gfa -V 2 -w map --prefix {DATASET_DIR}/sequence_map"
    subprocess.run(vg_command.split())
    
    # Convert xg index to vg for extracting the paths later
    #f = open(f"{DATASET_DIR}/sequence_map.vg", 'w')
    #vg_command = f"{vg_path} convert {DATASET_DIR}/sequence_map.xg -p"
    #subprocess.run(vg_command.split(), stdout=f)
    #f.close()
    print("Done")
