#!/bin/bash
#SBATCH --job-name=mbev024
#SBATCH --partition=largemem
#SBATCH --time=4-00:00:00
#SBATCH --mem=2500G
#SBATCH --qos=generous-compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=79
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

WORKDIR=$(pwd)
PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
query_genome="MB"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
#mkdir "$PATH_WGS"
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps/
ref_genome=EV_mazia
cpus=79
$APPSDIR/gfastats/build/bin/gfastats -f bananaB.subreads.fastq.gz -j $cpus -t > bananaB.subreads.txt

