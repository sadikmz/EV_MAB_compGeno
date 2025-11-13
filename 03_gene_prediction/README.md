### 1. Gene prediction 

- [BRAKER2](protein_coding_genes_prediction/braker.sh) with GeneMark-ETP and AUGUSTUS to predict *ab initio* gene models supported by protein homology and/or RNA-seq evidence. 

- RNA-seq aligned was generaed with [hisat2 (v2.1.0)](protein_coding_genes_prediction/hisat2.sh) using raw RNA-seq data of landrace Mazia. 

- RNAseq transcripts were assembeld with [Trinity (v2.1.1)](protein_coding_genes_prediction/trinity.sh)  further process to build splice-aware and unique transcript assemblies using [PASA pipeline (v2.4.1)](protein_coding_genes_prediction/pasa.sh)) 

- Final gene prediction was performed iteratively with [MAKER](rotein_coding_genes_prediction/trinity.sh) iteratively to train SNAP gene models and integrate it with protein homology hints and evidence from assembled RNA-seq transcripts.


### 2. Extract predicted gene models, filter gene models with selected AED value and assess BUSCO scores - [collect_filter_maker_annotation.sh](collect_filter_maker_annotation.sh) 

### 3. Functional annotation of predicted proteins - [functional_annotation.sh](functional_annotation.sh)

- Remove poor quality / TE-related gene models from final gene model set: [TEsorter_genemodels](TEsorter_genemodels.sh), [parse TE associated gene models](execlude_TE_inserted_gene_models.R) and [plot](parse_plot_TEsorter_out.R)