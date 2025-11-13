#!/bin/bash
#SBATCH --job-name=faniMA
#SBATCH --partition=standard
#SBATCH --time=14-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# data 
DATADIR="/home/data/salixomics/sadik/evmab_dir"
threads=32
ml FastANI

# run 

fastANI \
-r EMZ.fna \
-q EBD.fna \
-t $threads \
-o EBD_EMZ.fastani.out


