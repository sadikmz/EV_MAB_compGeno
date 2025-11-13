# extract glycerol-3-phosphate gene 

interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  filter(clust == "MA_specific" |
           clust == "MB_specific" |
           clust == "EV_mazia_specific" |
           clust == "EV_bedadeti_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      select(-GO_terms)
  ) %>%
  select(seqid, start, end, gene_id, clust, lastz_cov, GO_terms, description) %>%
  distinct() %>%
  # filter(#clust == "MA_specific" 
  #   #clust == "MB_specific" 
  #   clust == "EV_mazia_specific" 
  #     # clust == "EV_bedadeti_specific"
  # ) %>% 
  rename(seq.name = gene_id) %>%
  filter(#description=="carbohydrate metabolic process" |
    str_detect(description, "glycerol-3-phosphate") 
    # str_detect(description, "glycolytic process") 
    # str_detect(description, "defense")
    # str_detect(description, "glutathione metabolic process")
  ) %>%
  distinct(seqid,start,end) %>%
  write.table(paste0(result_dir, "EV_bedadeti_glycerol_3_phasphase_dehydrogenase.cds.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")


interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  filter(clust == "MA_specific" |
           clust == "MB_specific" |
           clust == "EV_mazia_specific" |
           clust == "EV_bedadeti_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      select(-GO_terms)
  ) %>%
  select(seqid, start, end, gene_id, clust, lastz_cov, GO_terms, description) %>%
  distinct() %>%
  # filter(#clust == "MA_specific" 
  #   #clust == "MB_specific" 
  #   clust == "EV_mazia_specific" 
  #     # clust == "EV_bedadeti_specific"
  # ) %>% 
  rename(seq.name = gene_id) %>%
  filter(#description=="carbohydrate metabolic process" |
    str_detect(description, "glycerol-3-phosphate") 
    # str_detect(description, "glycolytic process") 
    # str_detect(description, "defense")
    # str_detect(description, "glutathione metabolic process")
  ) %>%
  distinct() %>%
  
left_join(
  EV_MAB_prot_fasta %>%
    filter(str_detect(seq.name,"EVBD")) %>%
    distinct(seq.name) %>%
    separate(seq.name, into = c("seq.name", "isomer"), sep = "-") 
) %>% 
  mutate(seq_name = str_c(seq.name, isomer, sep = "-")) %>%
  distinct(seq_name) %>%
  rename(seq.name = seq_name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EVBD")) %>%
      select(-genome)
  ) %>%
  # rename(seq.name = seq_name) %>%
  # mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  # filter(seq.text == "NA")
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct()  %>%
  write.table(paste0(result_dir, "EV_bedadeti_glycerol_3_phasphase_dehydrogenase.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)
