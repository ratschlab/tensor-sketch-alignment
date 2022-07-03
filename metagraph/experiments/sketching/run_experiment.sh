rm *.png # cleanup

# Parameters
metagraph_path="/home/alex/metagraph/metagraph/build/metagraph"
output_path="./data"
max_k=90
graph_seq_len=100000000

num_query_seqs=1000
mutation_rate=10
min_path_size=200
max_path_size=201

parallel=$(getconf _NPROCESSORS_ONLN)
batch_size=1000
################

# Generate the dataset
gen_command="python generate_dataset.py --metagraph-path $metagraph_path --graph-seq-len $graph_seq_len --max-k $max_k"
eval $gen_command

# Sketch seeder params
embed_dim=64
minimizer_window=100
n_times_sketch=1
seeder="sketch"
######################

sketch_command="python seed_recall_on_k.py --minimizer-window $minimizer_window --output-path $output_path --embed-dim $embed_dim --n-times-sketch $n_times_sketch --mutation_rate $mutation_rate --num_query_seqs $num_query_seqs --parallel $parallel --batch-size $batch_size --seeder $seeder --metagraph-path $metagraph_path --max-k $max_k --max-path-size $max_path_size --min-path-size $min_path_size"

echo "[INFO] Launching the experiment on SKETCH"
eval $sketch_command

# Default seeder params
seeder="default"
#######################

default_command="python seed_recall_on_k.py --minimizer-window $minimizer_window --output-path $output_path --embed-dim $embed_dim --n-times-sketch $n_times_sketch --mutation_rate $mutation_rate --num_query_seqs $num_query_seqs --parallel $parallel --batch-size $batch_size --seeder $seeder --metagraph-path $metagraph_path --max-k $max_k --max-path-size $max_path_size --min-path-size $min_path_size"

echo "[INFO] Launching the experiment on DEFAULT"
eval $default_command

# Do overlap plot
runs=($(ls runs | sort | tail -n 2))
sketch_run="./runs/${runs[0]}/points.json"
default_run="./runs/${runs[1]}/points.json"

overlap_command="python overlap.py --max-k $max_k --sketch-json $sketch_run --exact-json $default_run"

echo $sketch_run
echo $default_run
echo "[INFO] Saving overlap plot"

eval $overlap_command

echo "Finished experiment"
