  EV_MAB_custered_genes.v5.1 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      distinct(gene_id, GO_terms)
  ) %>%
  # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
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
    write.table(paste0(result_dir, "MA_specific.CDS.bed"),
                col.names = F,
                row.names = F,
                sep = "\t",
                quote = F)
  
  
  ## MB
  
  EV_MAB_custered_genes.v5.1 %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        distinct(gene_id, GO_terms)
    ) %>%
    # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        select(-GO_terms)
    ) %>%
    distinct() %>%
    
    filter(clust == "MB_specific" 
           # clust == "MB_specific" 
           # clust == "EV_mazia_specific" 
           # clust == "EV_bedadeti_specific"
    ) %>% 
    # filter(str_detect(description,"systemic acquired resistance") |
    #          str_detect(description, "defense response")) %>%
    distinct(seqid, start, end, lastz_cov) %>%
    # nrow()
    # head()
    write.table(paste0(result_dir, "MB_specific.CDS.bed"),
                col.names = F,
                row.names = F,
                sep = "\t",
                quote = F)
  
  # EV mazia
  
  
  EV_MAB_custered_genes.v5.1 %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        distinct(gene_id, GO_terms)
    ) %>%
    # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        select(-GO_terms)
    ) %>%
    distinct() %>%
    
    filter(#clust == "MB_specific" 
           # clust == "MB_specific" 
           clust == "EV_mazia_specific"
           # clust == "EV_bedadeti_specific"
    ) %>% 
    # filter(str_detect(description,"systemic acquired resistance") |
    #          str_detect(description, "defense response")) %>%
    distinct(seqid, start, end, lastz_cov) %>%
    # nrow()
    # head()
    write.table(paste0(result_dir, "EV_maazia_specific.CDS.bed"),
                col.names = F,
                row.names = F,
                sep = "\t",
                quote = F)
  
  # EV bedadeti
  
  EV_MAB_custered_genes.v5.1 %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        distinct(gene_id, GO_terms)
    ) %>%
    # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        select(-GO_terms)
    ) %>%
    distinct() %>%
    
    filter(#clust == "MB_specific" 
      # clust == "MB_specific" 
      # clust == "EV_mazia_specific"
      clust == "EV_bedadeti_specific"
    ) %>% 
    # filter(str_detect(description,"systemic acquired resistance") |
    #          str_detect(description, "defense response")) %>%
    distinct(seqid, start, end, lastz_cov) %>%
    # nrow()
    # head()
    write.table(paste0(result_dir, "EV_bedadeti_specific.CDS.bed"),
                col.names = F,
                row.names = F,
                sep = "\t",
                quote = F)
  
  
  EV_MAB_custered_genes.v5.1 %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        distinct(gene_id, GO_terms)
    ) %>%
    # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        select(-GO_terms)
    ) %>%
    distinct() %>%
    
    filter(#clust == "MB_specific" 
      # clust == "MB_specific" 
      clust == "EV_mazia_specific"
      # clust == "EV_bedadeti_specific"
    ) %>% 
    # filter(str_detect(description,"systemic acquired resistance") |
    #          str_detect(description, "defense response")) %>%
    distinct(seqid, start, end, lastz_cov, gene_id) %>%
    # filter(lastz_cov == 0) %>%
    # distinct(seqid) %>%
    group_by(seqid) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    filter(count > 5) %>%
    arrange(desc(count)) %>%
    as.data.frame()

  
  
  ## Bedadeti
  
  EV_MAB_custered_genes.v5.1 %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>%
        distinct(gene_id, GO_terms)
    ) %>%
    # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms #%>%
        # select(-GO_terms)
    ) %>%
    distinct() %>%
    
    filter(#clust == "MB_specific" 
      # clust == "MB_specific" 
      # clust == "EV_mazia_specific"
      clust == "EV_bedadeti_specific"
    ) %>% 
    # filter(str_detect(description,"systemic acquired resistance") |
    #          str_detect(description, "defense response")) %>%
    distinct(seqid, start, end, lastz_cov, gene_id) %>%
    group_by(seqid) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    filter(count > 1) %>%
    arrange(desc(count)) %>%
    as.data.frame() %>%
    # head(n=50)
    # joining gene id and cds region 
    left_join(
      EV_MAB_custered_genes.v5.1 %>%
        left_join(
          EV_MAB_reads_coverage_geneID_GO_terms %>%
            distinct(gene_id, GO_terms)
        ) %>%
        # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
        left_join(
          EV_MAB_reads_coverage_geneID_GO_terms #%>%
          # select(-GO_terms)
        ) %>%
        distinct() %>%
        
        filter(#clust == "MB_specific" 
          # clust == "MB_specific" 
          clust == "EV_bedadeti_specific"
          # clust == "EV_bedadeti_specific"
        ) 
    ) %>%
    distinct(seqid, start, end, lastz_cov, gene_id, GO_terms) %>%
    filter(GO_terms != "NA") %>%
    left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T) %>%
    filter(aspects == "biological_process") %>%
    filter(!str_detect(description,"obsolete"))) %>%
    mutate(description = str_replace_na(description, "NA")) %>%
    filter(description != "NA")
    filter(description == "galactose metabolic process") #%>%
    distinct(description)
    head()
             
  
