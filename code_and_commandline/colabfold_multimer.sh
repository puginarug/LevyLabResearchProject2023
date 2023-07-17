#!/bin/bash

module load cuda/11.1
sbatch --gres gpu:a100-1-10 --time=36:0:0 --mem=24G --wrap="/sci/labs/asafle/alexlevylab/icore-data/tools/local_collabfold/v1.4.0/colabfold_batch/bin/colabfold_batch --num-recycle 15 --stop-at-score 0.90 --model-type AlphaFold2-multimer-v2 --rank multimer /sci/labs/asafle/davzilevylab/TA_project/fasta_to_fold/$1.fasta out_$1"
