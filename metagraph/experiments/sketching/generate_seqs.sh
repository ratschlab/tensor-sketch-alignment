#!/bin/bash

for i in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	python generate_seqs.py --dataset-dir /cluster/work/grlab/ameterez/sketch_data/data_$i
done
