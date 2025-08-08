EV_MAB_custered_genes.v1 %>%
  distinct(clust)

EV_MAB_custered_genes.v1 %>%
left_join(GO_terms_global_genes_MABST) %>%
  mutate(GO_terms = str_replace_na(GO_terms,"NA")) %>%
  filter(GO_terms != "NA") %>%
  distinct(gene_id,clust) %>%
  group_by(clust) %>%
  summarise(genes_count = n()) %>%
  write.table(paste0("Z:/sadik/evmab_dir/","EV_MAB_enrichment.gene_count.txt"),
              col.names = T, 
              row.names = F,
              quote = F,
              sep = "\t")
  

GO_terms_global_genes
EV_MAB_reads_coverage_geneID_GO_terms_wide.v1 

clustProfiler_GOenrich_output_MAST %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  write.table(paste0("Z:/sadik/evmab_dir/","EV_MAB_enrichment.MAST.txt"),
              col.names = T, 
              row.names = F,
              quote = F,
              sep = "\t")


# write enrichment output to disk 

clustProfiler_GOenrich_output %>%
  dplyr::select(-Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>% 
  write.table(paste0(DATA_DIR,"EV_MAB_clusterProfiler_enrichment_out.txt"),
              col.names = T,
              row.names = F, 
              sep = '\t',
              quote = F)
  # head()


clustProfiler_GOenrich_output %>%
  dplyr::select(-Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  # left_join(go_odb_aspects_description) %>%
  head()
