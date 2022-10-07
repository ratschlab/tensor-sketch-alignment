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
graph_seq_len=$1
num_levels=0
mutation_rate=0.01

mkdir data_$graph_seq_len
# Generate the dataset
gen_command="python generate_dataset.py --dataset-dir data_$graph_seq_len --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --mutation-rate $mutation_rate --num-levels $num_levels"
echo $gen_command
eval $gen_command
