############### Parameters ###############
# rm -rf data/*
metagraph_path="/home/alex/metagraph/metagraph/build/metagraph"
graph_seq_len=10000
num_levels=0
mutation_rate=0.05

# Generate the dataset
gen_command="python generate_dataset.py --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --mutation-rate $mutation_rate --num-levels $num_levels"
eval $gen_command
