### 1. Gene prediction 

- [BRAKER2](protein_coding_genes_prediction/braker.sh) with GeneMark-ETP and AUGUSTUS to predict *ab initio* gene models supported by protein homology and/or RNA-seq evidence. 

- RNA-seq aligned was generaed with [hisat2 (v2.1.0)](protein_coding_genes_prediction/hisat2.sh) using raw RNA-seq data of landrace Mazia. 

- RNAseq transcripts were assembeld with [Trinity (v2.1.1)](protein_coding_genes_prediction/trinity.sh)  further process to build splice-aware and unique transcript assemblies using [PASA pipeline (v2.4.1)](protein_coding_genes_prediction/pasa.sh)) 

- Final gene prediction was performed iteratively with [MAKER](rotein_coding_genes_prediction/trinity.sh) iteratively to train SNAP gene models and integrate it with protein homology hints and evidence from assembled RNA-seq transcripts.
- Predicted gene models were extraced and filtered for selected high-quality gene models - [collect_filter_maker_annotation.sh](collect_filter_maker_annotation.sh) 

