
# MA genes shared with EV
EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa_acuminata")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>%
  # head()
  left_join(
    EV_MAB_custered_genes.v5.1 %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          distinct(gene_id, GO_terms)
      ) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          select(-GO_terms)
      ) %>%
      distinct() %>%
      
      filter(clust == "MA_specific" 
             # clust == "MB_specific" 
             # clust == "EV_mazia_specific" 
             # clust == "EV_bedadeti_specific"
      ) %>% 
      distinct(seqid, start, end, gene_id) %>%
      # count here is to show the position on the genes in the plot
      mutate(count = 1)) %>%
  mutate(count = str_replace_na(count, "10")) %>%
  filter(count == "10") %>%
  mutate(count = str_remove(count,"0")) %>%
  dplyr::rename(seqnames = seqid) %>%
  dplyr::distinct(seqnames,start,end,count) %>%
  # filter(str_detect(seqnames,"chr")) %>%
  # mutate(seqnames = str_replace(seqnames, "chr", "A")) %>%
  # head()
  write.table(paste0(result_dir,"MA_genes_sharedWith_EV.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")



# MB genes shared with EV
EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa_acuminata")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>%
  # head()
  left_join(
    EV_MAB_custered_genes.v5.1 %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          distinct(gene_id, GO_terms)
      ) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          select(-GO_terms)
      ) %>%
      distinct() %>%
      
      filter(clust == "MA_specific" 
             # clust == "MB_specific"
             # clust == "EV_mazia_specific" 
             # clust == "EV_bedadeti_specific"
      ) %>% 
      distinct(seqid, start, end, gene_id) %>%
      # count here is to show the position on the genes in the plot
      mutate(count = 1)) %>%
  mutate(count = str_replace_na(count, "10")) %>%
  filter(count == "10") %>%
  mutate(count = str_remove(count,"0")) %>%
  dplyr::rename(seqnames = seqid) %>%
  dplyr::distinct(seqnames,start,end) %>%
  distinct() %>%
  # nrow()
  # nrow()
  # filter(!str_detect(seqnames,"chr")) %>%
  # nrow()
  write.table(paste0(result_dir,"MA_genes_sharedWith_EV.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")
  
  # write.table(paste0(result_dir,"MA_genes_sharedWith_EV.bed"),
  #             col.names = F,
  #             row.names = F,
  #             quote = F,
  #             sep = "\t")


EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa_balbisiana")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>%
  # head()
  left_join(
    EV_MAB_custered_genes.v5.1 %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          distinct(gene_id, GO_terms)
      ) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          select(-GO_terms)
      ) %>%
      distinct() %>%
      
      filter(#clust == "MA_specific" 
             clust == "MB_specific"
             # clust == "EV_mazia_specific" 
             # clust == "EV_bedadeti_specific"
      ) %>% 
      distinct(seqid, start, end, gene_id) %>%
      # count here is to show the position on the genes in the plot
      mutate(count = 1)) %>%
  mutate(count = str_replace_na(count, "10")) %>%
  filter(count == "10") %>%
  mutate(count = str_remove(count,"0")) %>%
  dplyr::rename(seqnames = seqid) %>%
  dplyr::distinct(seqnames,start,end, count) %>%
  distinct() %>%
  left_join(

read.delim(paste0(result_dir,"MB_sharedWith_EV_50kb_window.bed"), header=F) %>% 
  mutate(#V4 = case_when(str_detect(V4, "\\.") ~ "0",
          #              TRUE ~ "1"),
         #V4 = as.numeric(V4),
         # V3 = V3-1
         ) %>%
  # distinct(V1, V2, V3) %>%
  mutate (V8 = case_when(str_detect(V4, "\\w+") ~ "1",
                         TRUE ~ "0"),
          V8 = as.numeric(V8)) %>%
  # select(V1,V2,V3,V8) %>%
  # group_by(V1,V2,V3) %>%
  #   summarise(count = sum(V8)) %>%
  # ungroup() %>%
  # summarise(sum = sum(count))
  rename(seqnames =V4,
         start = V5,
         end = V6)) %>%
  mutate( count = str_replace_na(count, "NA")) %>%
  filter(count != "NA") %>% 
  distinct(seqnames,start,end, V1, V2, V3) %>%
  as.data.frame() %>%
  head(n=100)
  write.table(paste0(result_dir, "MB_50k_windows.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")


