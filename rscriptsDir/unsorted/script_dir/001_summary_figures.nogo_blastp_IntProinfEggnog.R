# EV_MAB specific genes that lack GO-terms annotation but showed blastp hits against ncbi nr, uniprot, and rich and arabidopsis proteom
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  head()


EV_MAB_specific_nogo_protinfer_interproscan_eggnog_MABST <- 
  EV_MAB_custered_genes.v2 %>% 
  filter(mapping_reads_source == "cross_species_reads") %>%
  # filter( query != "MS_LongReads",
  #           query != "MT_LongReads") %>%
  filter( frac_ovrerlap_syn < 0.245) %>%
  
  # filter(clust== "EV_mazia_specific") %>%
  filter(annot_intpro != "NA" |
           annot_eggnog != "NA" |
           annot_protenfer != "NA") %>%
  select(seqid, cds_start, cds_end, genome, gene_id, predicted_label, protinfer_description,eggNOG_OGs,Preferred_name) %>%
  mutate(#interproscan_description = str_replace_na(interproscan_description,"NA"),
    protinfer_description = str_replace_na(protinfer_description,"NA"),
    eggNOG_OGs = str_replace_na(eggNOG_OGs,"NA"),
    Preferred_name = str_replace_na(Preferred_name,"NA"))  %>%
  distinct(gene_id,predicted_label,protinfer_description) %>%
  filter(protinfer_description != "NA")