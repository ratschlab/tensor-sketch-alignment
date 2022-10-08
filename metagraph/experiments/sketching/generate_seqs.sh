#!/bin/bash

for i in 21 22
do
	python generate_seqs.py --dataset-dir /cluster/work/grlab/ameterez/sketch_data/data_$i
done
