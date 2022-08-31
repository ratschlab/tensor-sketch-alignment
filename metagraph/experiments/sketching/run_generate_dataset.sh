############### Parameters ###############
# rm -rf data/*
metagraph_path="/home/alex/metagraph/metagraph/build/metagraph"
max_k=81
graph_seq_len=100000
num_seqs=1
mutation_rate=0.05

# Generate the dataset
gen_command="python generate_dataset.py --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --max-k $max_k --mutation-rate $mutation_rate --num-seqs $num_seqs"
eval $gen_command
