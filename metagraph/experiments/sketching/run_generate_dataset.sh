############### Parameters ###############
# rm -rf data/*
#
#
#
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi
metagraph_path="/cluster/apps/biomed/grlab/ameterez/metagraph/metagraph/build/metagraph"
graph_seq_len=1000
num_levels=$1
mutation_rate=0.01

mkdir data_$num_levels
# Generate the dataset
gen_command="python generate_dataset.py --dataset-dir data_$num_levels --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --mutation-rate $mutation_rate --num-levels $num_levels"
echo $gen_command
eval $gen_command
