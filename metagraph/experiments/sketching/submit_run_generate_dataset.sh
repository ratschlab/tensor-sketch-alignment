#!/bin/bash
sbatch << EOT
#!/bin/bash

#SBATCH -o "run_$1.stdout"
#SBATCH -e "run_$1.stderr"
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=60G

source ~/.bashrc
enable_modules
conda activate base

mkdir /cluster/work/grlab/ameterez/sketch_data/data_$1
# Generate the dataset
python generate_dataset.py --dataset-dir /cluster/work/grlab/ameterez/sketch_data/data_$1 --metagraph-path /cluster/apps/biomed/grlab/ameterez/metagraph/metagraph/build/metagraph --graph-seq-len 1200 --mutation-rate 0.01 --num-levels $1 

exit 0
EOT
