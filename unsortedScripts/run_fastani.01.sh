#!/bin/bash
#SBATCH --job-name=faniMA
#SBATCH --partition=largemem
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# data 
DATADIR="/home/data/salixomics/sadik/evmab_dir"
threads=80
ml FastANI


# run 

fastANI \
-r MA.fna \
-q EBD.fna \
-t $threads \
-o EBD_MA.fastani.00.out

fastANI \
-r MB.fna \
-q EBD.fna \
-t $threads \
-o EBD_MB.fastani.00.out

