#!/bin/bash

module load cuda/11.1
sbatch --gres gpu:a100-1-10 --time=5:0:0 --mem=10G --wrap="/sci/labs/asafle/alexlevylab/icore-data/tools/local_collabfold/v1.4.0/colabfold_batch/bin/colabfold_batch /sci/labs/asafle/davzilevylab/TA_project/fasta_to_fold/$1.fasta out_$1"
