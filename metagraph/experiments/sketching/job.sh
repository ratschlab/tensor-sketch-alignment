
#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -o "runs/run_method_$1_mutation_$2_data_$3_suffix_$4.stdout"
#SBATCH -e "runs/run_method_$1_mutation_$2_data_$3_suffix_$4.stderr"
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=80G

source ~/.bashrc
enable_modules
conda activate base

python run_one_baseline.py --method $1 --mutation $2 --dataset-dir /cluster/work/grlab/ameterez/sketch_data/data_$3 --suffix $4

exit 0
EOT
