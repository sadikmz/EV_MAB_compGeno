#!/bin/bash
#SBATCH --job-name=0041
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


module load InterProScan_data/5.55-88.0-foss-2021a


DATADIR="/home/data/salixomics/sadik/evmab_dir"
APPSDIR="/home/data/salixomics/sadik/apps/interproscan-5.69-101.0"
prefix=EV_MAB_nohit_genes
threads=80

cat EV_MAB_nogo_noblastp_noproteinfer.fa | sed 's/*//g' > EV_MAB_nogo_noblastp_noproteinfer.01.fa

"$APPSDIR"/interproscan.sh \
-cpu $threads \
-pa \
-i "$DATADIR"/EV_MAB_nogo_noblastp_noproteinfer.01.fa \
-f TSV,GFF3

