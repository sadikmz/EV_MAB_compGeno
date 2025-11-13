# Write to disk GO terms associated with biological process for each speceis

interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "MA_specific") %>% 
  dplyr::filter(description != "NA") %>%
  dplyr::rename(MA_gene_list = gene_list) %>%
  select(category, description, Count,MA_gene_list) %>%
  mutate(category = str_remove(category,"_specific")) %>%
  dplyr::rename(#cat = category,
    MA = Count) %>%
  full_join(
    
    interproscan_GO_terms.v1 %>%
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "MB_specific") %>%
      dplyr::filter(description != "NA") %>%
      dplyr::rename(MB_gene_list = gene_list) %>%
      select(category, description, Count, MB_gene_list) %>%
      mutate(category = str_remove(category,"_specific")) %>%
      dplyr::rename(cat = category,
             MB = Count)) %>%
  dplyr::select(description,MA,MB, MA_gene_list, MB_gene_list) %>%
  mutate(MA = str_replace_na(MA, "0"),
         MB = str_replace_na(MB, "0"),
         MA = as.numeric(MA),
         MB = as.numeric(MB)) %>%
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T) %>%
      filter(aspects == "biological_process") %>%
      filter(!str_detect(description,"obsolete")) 
  ) %>% 
  dplyr::select(GO_terms,description,MA,MB,MA_gene_list,MB_gene_list) %>%
  # distinct(GO_terms,description,MA,MB) %>%
  # 
  # filter(str_detect(description,"phospho") |
  #          str_detect(description, "signal transduction") |
  #          str_detect(description, "regulation")) 

  write.table(paste0(result_dir, "MAB_go_terms_biological_process.txt"),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)


## For EV mazia 

interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "EV_mazia_specific") %>% 
  dplyr::filter(description != "NA") %>%
  dplyr::rename(EV_mazia_gene_list = gene_list) %>%
  select(category, description, Count,EV_mazia_gene_list) %>%
  mutate(category = str_remove(category,"_specific")) %>% 
  dplyr::rename(#cat = category,
    EV_mazia = Count) %>%
  full_join(
    
    interproscan_GO_terms.v1 %>%
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "EV_bedadeti_specific") %>%
      dplyr::filter(description != "NA") %>%
      dplyr::rename(EV_bedadeti_gene_list = gene_list) %>%
      select(category, description, Count, EV_bedadeti_gene_list) %>%
      mutate(category = str_remove(category,"_specific")) %>%
      dplyr::rename(cat = category,
             EV_bedadeti = Count)) %>%
  select(description,EV_mazia,EV_bedadeti, EV_mazia_gene_list, EV_bedadeti_gene_list) %>%
  mutate(EV_mazia = str_replace_na(EV_mazia, "0"),
         EV_bedadeti = str_replace_na(EV_bedadeti, "0"),
         EV_mazia = as.numeric(EV_mazia),
         EV_bedadeti = as.numeric(EV_bedadeti)) %>%
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T) %>%
      filter(aspects == "biological_process") %>%
      filter(!str_detect(description,"obsolete")) 
  ) %>% 
  select(GO_terms,description,EV_mazia,EV_bedadeti,EV_mazia_gene_list,EV_bedadeti_gene_list) %>%
  distinct() %>%
  # distinct(GO_terms,description,EV_mazia,EV_bedadeti) %>%
  # filter(str_detect(description,"phospho") |
  #          str_detect(description, "signal") |
  #          str_detect(description, "regulation")) 

  write.table(paste0(result_dir, "EV_go_terms_biological_process.txt"),
              col.names = T,
              row.names = F,
              sep = "\t",
              quote = F)

# GO_terms CDS

# interproscan_GO_terms.v1 <-
  interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
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
  filter(#clust == "MA_specific" 
    clust == "MB_specific" 
    # clust == "EV_mazia_specific" 
    # clust == "EV_bedadeti_specific"
  ) %>% 
  distinct(seqid, start, end, description, gene_id) %>%
  write.table(paste0(result_dir, "MB_go_terms_biological_process.CDS.bed"),
              col.names = F,
              row.names = F,
              sep = "\t",
              quote = F)

# interproscan_GO_terms.v1 <-
interproscan_GO_terms %>% 
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
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
  distinct() %>%
  filter(clust == "MA_specific" 
         # clust == "MB_specific" 
         # clust == "EV_mazia_specific" 
         # clust == "EV_bedadeti_specific"
  ) %>% 
  # filter(str_detect(description,"systemic acquired resistance") |
  #          str_detect(description, "defense response")) %>%
  distinct(seqid, start, end, lastz_cov) %>%
  # head()
write.table(paste0(result_dir, "MA_go_terms_biological_process.CDS.bed"),
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)

interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
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
  distinct() %>%
  filter(clust == "MA_specific" 
    # clust == "MB_specific" 
    # clust == "EV_mazia_specific" 
    # clust == "EV_bedadeti_specific"
  ) %>% 
  filter(str_detect(description,"systemic acquired resistance") |
           str_detect(description, "defense response")) %>%
  distinct(seqid, start, end, description) %>%
  write.table(paste0(result_dir, "MA_go_terms_BP_disease_resistance.CDS.bed"),
              col.names = F,
              row.names = F,
              sep = "\t",
              quote = F)


interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
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
  distinct() %>%
  filter(#clust == "MA_specific" 
    clust == "MB_specific" 
    # clust == "EV_mazia_specific" 
    # clust == "EV_bedadeti_specific"
  ) %>% 
  filter(str_detect(description,"systemic acquired resistance") |
           str_detect(description, "defense response")) %>%
  select(seqid, start, end, lastz_cov ) %>%
  write.table(paste0(result_dir, "MB_go_terms_BP_disease_resistance.CDS.bed"),
              col.names = F,
              row.names = F,
              sep = "\t",
              quote = F)


# protein coding genes 
## EV mazia
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  filter(#clust == "MA_specific" 
    #clust == "MB_specific" 
    clust == "EV_mazia_specific"
    # clust == "EV_bedadeti_specific"
  ) %>% 
  rename(seq.name = gene_id) %>%
  distinct(seq.name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EVMZ")) %>%
      distinct(seq.name) %>%
      separate(seq.name, into = c("seq.name", "isomer"), sep = "-") 
  ) %>% 
  mutate(seq_name = str_c(seq.name, isomer, sep = "-")) %>%
  distinct(seq_name) %>%
  rename(seq.name = seq_name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EVMZ")) %>%
      select(-genome)
  ) %>%
  # rename(seq.name = seq_name) %>%
  # mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  # filter(seq.text == "NA")
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>% 
  # head()
  write.table(paste0(result_dir, "EV_mazia_go_terms_biological_process.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)


# cds

interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
  filter(clust == "MA_specific" |
           clust == "MB_specific" |
           clust == "EV_mazia_specific" |
           clust == "EV_bedadeti_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      select(-GO_terms)
  ) %>%
  select(seqid, start, end, lastz_cov, gene_id, clust) %>%
  distinct() %>%
  filter(#clust == "MA_specific" 
    #clust == "MB_specific" 
    clust == "EV_mazia_specific"
    # clust == "EV_bedadeti_specific"
  ) %>% 
  distinct(seqid, start, end, lastz_cov, gene_id) %>%
  arrange(seqid) %>%
  head()
  write.table(paste0(result_dir, "EV_mazia_go_terms_biological_process.cds.bed"),
              col.names = F,
              row.names = F,
              sep = "\t",
              quote = F)

## EV_bedadeti

interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  filter(#clust == "MA_specific" 
    #clust == "MB_specific" 
    # clust == "EV_mazia_specific"
    clust == "EV_bedadeti_specific"
  ) %>% 
  rename(seq.name = gene_id) %>%
  distinct(seq.name) %>%
  # head()
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
  distinct() %>% 
  # tail()
  write.table(paste0(result_dir, "EV_bedadeti_go_terms_biological_process.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)




# bed file
  
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>% 
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  filter(#clust == "MA_specific" 
    #clust == "MB_specific" 
    # clust == "EV_mazia_specific"
    clust == "EV_bedadeti_specific"
  ) %>% 
  # dplyr::rename(seq.name = gene_id) %>% 
  distinct(seqid, start, end, lastz_cov) %>% 
  write.table(paste0(result_dir, "EV_bedadeti_go_terms_biological_process.cds.bed"),
              col.names = F,
              row.names = F,
              sep = "\t",
              quote = F)


## MA
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  filter(clust == "MA_specific" 
    # clust == "MB_specific" 
    # clust == "EV_mazia_specific" 
    # clust == "EV_bedadeti_specific"
  ) %>% 
  rename(seq.name = gene_id) %>%
  distinct(seq.name) %>%
  left_join(EV_MAB_prot_fasta) %>% 
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>%
  write.table(paste0(result_dir, "MA_go_terms_biological_process.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)
  
## MB
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  filter(#clust == "MA_specific" 
    clust == "MB_specific" 
    # clust == "EV_mazia_specific" 
    # clust == "EV_bedadeti_specific"
  ) %>% 
  rename(seq.name = gene_id) %>%
  distinct(seq.name) %>%
  left_join(EV_MAB_prot_fasta) %>% 
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>%
  # head(n=10) %>%
  write.table(paste0(result_dir, "MB_go_terms_biological_process.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)


## Specific go terms genes 
## EV mazia
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  # select(gene_id,description,clust) %>% 
  # # distinct() 
  # group_by(clust,description) %>%
  # summarise(total_genes = n(),
  #           gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  # as.data.frame() %>%
  # arrange(desc(total_genes)) %>%
  # as.data.frame() %>% 
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
  #   clust == "EV_mazia_specific"|
  #   clust == "EV_bedadeti_specific"
  # ) %>% 
  rename(seq.name = gene_id) %>%
  # distinct(description,GO_terms)
  filter(str_detect(description, "carbohydrate metabolic process") )%>%
  # distinct(seq.name) %>%
  # head() %>%
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EV")) %>%
      distinct(seq.name) %>%
      separate(seq.name, into = c("seq.name", "isomer"), sep = "-") 
  ) %>% 
  mutate(
    seq_name = case_when(str_detect(seq.name,"EV") ~ str_c(seq.name, isomer, sep = "-"),
                         TRUE ~ seq.name)
    
    
    # seq_name = str_c(seq.name, isomer, sep = "-")
    
    ) %>%
  distinct(seq_name) %>%
  rename(seq.name = seq_name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      # filter(str_detect(seq.name,"EV")) %>%
      select(-genome)
  ) %>%
  # rename(seq.name = seq_name) %>%
  # mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  # filter(seq.text == "NA")
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>% 
  write.table(paste0(result_dir, "EV_MAB_carbo_metabolic_process.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)




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
  filter(#clust == "MA_specific" 
    #clust == "MB_specific" 
    clust == "EV_mazia_specific" |
    clust == "EV_bedadeti_specific"
  ) %>% 
  rename(seq.name = gene_id) %>%
  # distinct(description,GO_terms)
  filter(str_detect(GO_terms, "GO:0015074") |
           str_detect(GO_terms, "GO:0006313"))%>%
  # distinct(seq.name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EV")) %>%
      distinct(seq.name) %>%
      separate(seq.name, into = c("seq.name", "isomer"), sep = "-") 
  ) %>% 
  mutate(seq_name = str_c(seq.name, isomer, sep = "-")) %>%
  distinct(seq_name) %>%
  rename(seq.name = seq_name) %>%
  # head()
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EV")) %>%
      select(-genome)
  ) %>%
  # rename(seq.name = seq_name) %>%
  # mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  # filter(seq.text == "NA")
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>% 
  # head()
  write.table(paste0(result_dir, "EV_DNA_integration_transposition.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)
