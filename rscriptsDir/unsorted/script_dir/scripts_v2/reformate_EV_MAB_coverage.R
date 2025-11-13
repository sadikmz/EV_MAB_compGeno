



  #panEV reads against MA 
panEV_MA_frac_ovl_syn <-
  EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
  filter(query=="panEV_reads") %>%
  filter(ref_genome == "musa_acuminata") %>%
  filter(str_detect(seqid,"chr")) %>%
  select(seqid,start,end,frac_ovrerlap_syn) %>%
  mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) %>%
  distinct() 

MA_gene_bed <-
  EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
  filter(query=="panEV_reads") %>%
  filter(ref_genome == "musa_acuminata") %>%
  filter(str_detect(seqid,"chr")) %>%
  select(seqid,start,end) %>%
  distinct() 

  MA_gene_complement <-
    read.delim(paste0(DATA_DIR,"Musa_acuminata.gene.complement.bed"), header = F) #%>%
      # mutate(frac_ovrerlap_syn = 0.2)

  
  panEV_MB_frac_ovl_syn <-
    EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
    filter(query=="panEV_reads") %>%
    filter(ref_genome == "musa_balbisiana") %>%
    filter(str_detect(seqid,"chr")) %>%
    select(seqid,start,end,frac_ovrerlap_syn) %>%
    mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) %>%
    distinct() 
  
  MB_gene_bed <-
    EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
    filter(query=="panEV_reads") %>%
    filter(ref_genome == "musa_balbisiana") %>%
    filter(str_detect(seqid,"chr")) %>%
    select(seqid,start,end) %>%
    distinct() 
  
  MB_gene_bed %>%
    head()
  
  MB_gene_complement <-
    read.delim(paste0(DATA_DIR,"Musa_balbisiana.gene.complement.bed"), header = F) 
  
  MB_gene_complement %>%
    head()
  
  

MA_specific_genes <-
    EV_MAB_custered_genes.v1 %>%
        distinct() %>% 
        filter(str_detect(clust,"MA_specific")) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
          filter(query=="panEV_reads") %>%
          filter(ref_genome == "musa_acuminata") %>%
          filter(str_detect(seqid,"chr")) %>%
          select(seqid,start,end,gene_id,frac_ovrerlap_syn) %>%
          mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) 
      ) %>%
      distinct() %>%
      select(seqid,start,end,frac_ovrerlap_syn)

    
    # head()
    # 
    # select(seqid, cds_start,cds_end, clust ) %>%
    # head()
    # select(seqid, cds_start, cds_end, genome, gene_id, reads_genome, query, mapping_reads_source, fraction_overlap, syn_cov, frac_ovrerlap_syn, annot_no_blastp,clust )
    # 
    # 
  
EV_MAB_custered_genes.v1 %>%
  distinct() %>% 
  filter(str_detect(clust,"MA_specific")) %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
      filter(query=="panEV_reads") %>%
      filter(ref_genome == "musa_acuminata") %>%
      filter(str_detect(seqid,"chr")) %>%
      select(seqid,start,end,gene_id,frac_ovrerlap_syn) %>%
      mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) 
  ) %>%
  distinct() %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>%
  write.table(paste0(DATA_DIR,"Musa_acuminata.specific_genes_all.bed"), 
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t') 

  head()

    MB_specific_genes <-
      EV_MAB_custered_genes.v1 %>%
      distinct() %>% 
      filter(str_detect(clust,"MB_specific")) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
          filter(query=="panEV_reads") %>%
          filter(ref_genome == "musa_balbisiana") %>%
          filter(str_detect(seqid,"chr")) %>%
          select(seqid,start,end,gene_id,frac_ovrerlap_syn) %>%
          mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) 
      ) %>%
      distinct() %>%
      select(seqid,start,end,frac_ovrerlap_syn)
  
  
    EV_MAB_custered_genes.v1 %>%
      distinct() %>% 
      filter(str_detect(clust,"MB_specific")) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1 %>%
          filter(query=="panEV_reads") %>%
          filter(ref_genome == "musa_balbisiana") %>%
          filter(str_detect(seqid,"chr")) %>%
          select(seqid,start,end,gene_id,frac_ovrerlap_syn) %>%
          mutate(frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) 
      ) %>%
      mutate(seqid = str_replace_na(seqid,"scaffolds")) %>%
      distinct() %>%
      filter(seqid != "scaffolds") %>%
      select(seqid,start,end,gene_id) %>%
      distinct() %>% 
      write.table(paste0(DATA_DIR,"Musa_balbisiana.specific_genes_all.bed"), 
                  col.names = F,
                  row.names = F,
                  quote = F,
                  sep = '\t')
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  