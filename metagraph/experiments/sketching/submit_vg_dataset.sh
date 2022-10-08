#!/bin/bash
sbatch << EOT
#!/bin/bash

#SBATCH --job-name=vg_indexing_$1
#SBATCH -o "run_vg_$1.stdout"
#SBATCH -e "run_vg_$1.stderr"
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=80G

source ~/.bashrc
enable_modules
conda activate base

mkdir /cluster/work/grlab/ameterez/sketch_data/data_$1
# Generate the dataset
python vg_generate_dataset.py --dataset-dir /cluster/work/grlab/ameterez/sketch_data/data_$1
exit 0
EOT
