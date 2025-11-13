#!/bin/bash
#SBATCH --job-name=np_nwc1059
#SBATCH --partition=largemem
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
module load RAxML/8.2.12-gompi-2020a-hybrid-avx2
threads=80

# /home/data/salixomics/sadik/evmab_dir/phylogene_dir
apps_dir="/home/data/salixomics/sadik/apps/trimal-1.5.0/source"
# NLR
clustalo \
--in  EV_MAB.novel_vs_zingibrales_plastome.01.fa \
--out EV_MAB_novel_genes.01.phy \
--outfmt=phylip \
--threads $threads \
--seqtype=Protein

#### Trim alignments 
$apps_dir/trimal -in EV_MAB_novel_genes.01.phy -out EV_MAB_novel_genes_gt_90.phy -phylip -gt 0.90 
#positions with gaps above 90% in aligned sequences were removed using 

$apps_dir/trimal -in EV_MAB_novel_genes.01.phy -out EV_MAB_novel_genes_gt_80.phy -phylip -gt 0.80 
# #positions with gaps above 20% in aligned sequences were removed using 

$apps_dir/trimal -in EV_MAB_novel_genes.01.phy -out EV_MAB_novel_genes_gt_50.phy -phylip -gt 0.50 
#positions with gaps above 50% in aligned sequences were removed using 

# run raxml
raxmlHPC-HYBRID-AVX2 -T $threads -s EV_MAB_novel_genes_gt_90.phy -n EV_MAB_novel_genes_90_RAxML -m PROTGAMMAAUTO -f a -# 100 -x 123456 -p 246810 
#raxmlHPC-HYBRID-AVX2 -T $threads -s EV_MAB_novel_genes_gt_80.phy -n EV_MAB_novel_genes_80_RAxML -m PROTGAMMAAUTO -f a -# 100 -x 123456 -p 246810 
#raxmlHPC-HYBRID-AVX2 -T $threads -s EV_MAB_novel_genes_gt_50.phy -n EV_MAB_novel_genes_50_RAxML -m PROTGAMMAAUTO -f a -# 100 -x 123456 -p 246810 


