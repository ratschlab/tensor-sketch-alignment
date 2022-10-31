# Aligning Distant Sequences to Graphs using Long Seed Sketches

Code for the paper:

> __Aligning Distant Sequences to Graphs using Long Seed Sketches__  
> Amir Joudaki*, Alexandru Meterez*, Harun Mustafa, Ragnar Groot Koerkamp, Andre Kahles, Gunnar Raetsch

> __bioRxiv 2022.10.26.513890; doi: https://doi.org/10.1101/2022.10.26.513890__

## Documentation

Instalation and commands on how to reproduce the experiments in the paper and use our software.

### Install from source
The code is built on top of [MetaGraph](https://metagraph.ethz.ch/). 
To install the package from source, follow the instructions from the [MetaGraph documentation](https://metagraph.ethz.ch/static/docs/installation.html#install-from-source).

Notice that during installation, after initializing the git submodules (between steps 4 and 5), run the commands:
```
cd metagraph/external-libraries/faiss
mkdir build
cd build
cmake ..
make -j
```
This will build FAISS, which is a required dependency of the package.

### Usage

- Build a De Bruijn graph
```
./metagraph build 
            -k 80 
            $FASTA
```

- Generate random walk in the graph
```
./metagraph seqgen 
            --output-path $DATASET_DIR
            --mutation-rate $MR
            --min-path-size $MINPATH
            --max-path-size $MAXPATH
            --num-query-seqs $NUM_SEQS
            -i $DBG_PATH
            --experiment
```

- Run the Tensor Sketch seeder and alignment:
```
./metagraph align
            --seeder sketch
            --embed-dim $EMBED_DIM
            --n-times-sketch 1
            --parallel $NUM_THREADS
            --batch-size 1000
            --num-neighbours $KNN
            --align-end-bonus 0
            --load-index 0
            -i $DBG_PATH
            $FASTA
```

For further tunable parameters on the alignment, consult [the MetaGraph documentation](https://github.com/ratschlab/tensor-sketch-alignment).

### Reproduce Experiments
Detailed explanations will be added soon. The files that have been run to create the plots are found in `metagraph/experiments/sketching`.

## Developer Notes
Metagraph and the Tensor Sketch Aligner is distributed under the GPLv3 License (see LICENSE). Please find further information in the AUTHORS and COPYRIGHTS files.