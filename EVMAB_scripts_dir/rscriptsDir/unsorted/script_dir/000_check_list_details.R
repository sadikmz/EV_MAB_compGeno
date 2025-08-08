# from 00_ev_mab_mapping_cov

## Reads coverage    
EV_MAB_reads_coverage 

## Get regions in Enset genome syntenic to MA and MB that were inferrred from lastz and Mummer alignment 

## lastz and nucmer alignment coverage
EV_MAB_lastz_nucmer_coding_seq

## NGS reads and nucmer+lastz combined 
EV_MAB_reads_coverage_no_lastz_nucmer_alignment

# from 02_reads_coverage_go_terms
EV_MAB_genes
EV_MAB_reads_coverage_geneID
EV_MAB_reads_coverage_geneID_GO_terms
EV_MAB_reads_coverage_geneID_synteny %>% head()
EV_MAB_reads_coverage_geneID_GO_terms_synteny %>% head()

## Filtered only for genic regions, used 0.245 cutoff and added longreds of Musa species 
EV_MAB_reads_coverage_no_lastz_nucmer_alignment.v3
EV_MAB_reads_coverage_no_lastz_nucmer_alignment
EV_MAB_reads_coverage_geneID_synteny.v1


# From 04_panEV_MAB_specific_gene_enrichment_allMapped_non_syntenic.01

EV_MAB_reads_coverage_geneID_GO_terms_wide.v1  %>% head()
EV_MAB_reads_coverage_geneID_GO_terms_long 
EV_MAB_custered_genes.v2 %>% head
  distinct(clust)


# From 05_EV_MAB_cluterprofiler_enrichment.synteny

EV_MAB_protinfer
EV_MAB_protinfer_retrieved
eggnog_MAB_GO_terms
eggnog_MAB_retrieved
EV_MAB_eggnog_retrieved

EV_MAB_specific_genes_nogoterms
GO_terms_global_genes
clustProfiler_GOenrich_output



GO_terms_global_genes %>%
  mutate(genome = case_when(str_detect(gene_id,"EVMZ") ~ "EV_mazia",
                            str_detect(gene_id,"EVBD") ~ "EV_bedadeti",
                            str_detect(gene_id,"Mac") ~ "MA",
                            str_detect(gene_id,"Mb") ~"MB"
  )) %>% 
  left_join(go_odb_aspects_description) %>% 
  select(gene_id,genome,aspects) %>%
  distinct() %>% 
  group_by(genome,aspects) %>%
  summarise(count = n())
# From 05_EV_MABST_cluterprofiler_enrichment.synteny

EV_MAB_eggnog_retrieved_MABST
EV_MAB_specific_genes_nogoterms_MABST

GO_terms_global_genes_MABST %>%
clustProfiler_GOenrich_output_MAST
bedadeti_enrich_MABST_plot
bedadeti_enrich_MABST_plot.v1
mazia_enrich_MABST
mazia_enrich_MABST_plot
mazia_enrich_MABST_plot.v1
MA_enrich_MABST
MA_enrich_MABST_plot
MA_enrich_MABST_plot.v1
MB_enrich_MBAST
MB_enrich_MBAST_plot
MB_enrich_MBAST_plot.v1

# From 06_extract_nogo_seq

EV_MAB_prot_fasta
EV_MAB_specific_genes_nogoterms
GO_terms_EV_AED25_MAB_desc.v1

# From 07_EV_MAB_syn_nogo_genes_ncbi_uniprot and 08_EV_MABST_syn_nogo_genes_at_os

EV_MAB_nogo_noreads_ncbi_nr
EV_MAB_nogo_noreads_uniprot_nr
EV_MAB_nogo_noreads_at_os_nr
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits
EV_MAB_nogo_noreads_ncbi_uniprot_hits
EV_MAB_specific_genes_nogo_noblastp %>%
  head()

EV_MAB_specific_genes_nogo_noblastp_MAST %>%
  head()
EV_MAB_specific_nogo_protinfer_interproscan_eggnog_MABST
# From 11_interproscan_nogo_noblastp
EV_MAB_interprocan

# 13_EVMAB_novel_genes.R
EV_MAB_custered_genes.v2 
EV_MAB_custered_genes.v3
