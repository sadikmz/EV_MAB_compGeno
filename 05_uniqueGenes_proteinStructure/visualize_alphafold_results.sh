#!/bin/bash

module load Python/3.10.8-GCCcore-12.2.0
module load matplotlib

NAME=rbd
DIR=${NAME}/${NAME}/

python ./visualize_alphafold_results.py --input_dir $DIR --output_dir $DIR  ##--name $NAME
echo Done visualize_alphafold_results for $NAME 

