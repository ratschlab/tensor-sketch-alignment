#!/bin/bash

for dataset_dir in 500 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000
do
	for method in "map" "mpmap" "sketch" "default" "graphaligner"
	do
		for mutation in 0 5 10 15 20 25
		do
			bash job.sh $method $mutation $dataset_dir run_1
		done
	done
done
