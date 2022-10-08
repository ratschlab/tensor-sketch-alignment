#!/bin/bash

for method in "sketch" #"default" "graphaligner"
do
	for dataset_dir in {19..21}
	do
		for mutation in 0 5 10 15 20 25
		do
			bash job.sh $method $mutation $dataset_dir run_2
		done
	done
done
