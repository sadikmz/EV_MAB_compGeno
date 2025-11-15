#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --output=%j.%x.out
#SBATCH --error=%j.%x.err
#SBATCH --partition=standard
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --qos=medium
###SBATCH --mail-type=END,FAIL
###SBATCH --mail-user=an.other@rothamsted.ac.uk
set -eu 
echo hostname=$(hostname)
echo pwd=$(pwd)
echo USER=$USER
date +"%d-%b-%y %T"
#====================
#
##	based on https://elearning.vib.be/courses/alphafold/lessons/alphafold-on-the-hpc/topic/setting-up-the-job-script/

module load AlphaFold
export ALPHAFOLD_DATA_DIR=/home/data/bioinf_resources/Alphafold/20231026

PROTEIN=rbd.faa
NAME=$(basename $PROTEIN .faa)
OUTDIR=./${NAME}
mkdir -p $OUTDIR
PRESET=monomer_ptm	#or multimer

alphafold --fasta_paths=${PROTEIN} --max_template_date=2999-01-01 --db_preset=full_dbs \
	--output_dir=${OUTDIR} --model_preset=${PRESET} 

echo "Done SLURM job: "$SLURM_JOBID
echo "Done running" $0
date +"%d-%b-%y %T"
exit

