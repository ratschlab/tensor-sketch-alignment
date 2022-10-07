#!/bin/bash

for method in "sketch" "default" "graphaligner"
do
	for dataset_dir in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
	do
		for mutation in 0 5 10 15 20 25
		do
			bash job.sh $method $mutation $dataset_dir run_1
		done
	done
done
