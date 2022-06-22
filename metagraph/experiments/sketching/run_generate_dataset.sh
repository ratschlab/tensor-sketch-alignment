############### Parameters ###############
metagraph_path="/Users/alex/metagraph/metagraph/build/metagraph"
max_k=30
graph_seq_len=100000

# Generate the dataset
gen_command="python generate_dataset.py --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --max-k $max_k"
eval $gen_command
