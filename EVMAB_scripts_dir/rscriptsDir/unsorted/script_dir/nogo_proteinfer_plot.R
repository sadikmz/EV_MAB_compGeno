# GO_terms retrived from proteinfer annotation 

EV_MAB_protinfer %>%
  rename(gene_id = sequence_name) %>%
  filter(!str_detect(predicted_label,"GO:")) %>%
  filter(confidence > 0.9 ) %>%
  separate(gene_id,into = c("gene_id","isoform"),sep = "\\-") %>%
  head()


EV_MAB_protinfer_retrieved

# EV_MAB specific genes that do not have GO term but showed blastp against protein databases 


EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  head()


#EV_MAB se

# EV_MBA gene clusters 
  
  
EV_MAB_custered_genes.v1 %>%
  filter(clust=="EV_bedadeti_specific" |
           clust=="EV_mazia_specific" |
           clust=="MA_specific"|
           clust=="MB_specific" ) %>% 
  group_by(clust) %>%
  summarise(count = n())
  head()
  

# EV_MAB specific genes who have Go-terms 

GO_terms_global_genes %>%
  head()
  

  # Final genes novel EV_MAB genes list
  
EV_MAB_custered_genes.v3 %>% 
  group_by(genome) %>%
  summarise(count = n())
  nrow()

  
# Extract EV_MAB specific genes that lack GO-terms and blastp-based homology 
  
  proteinfer_nogo_noblastp_00 <- 
  EV_MAB_custered_genes.v1 %>%
    filter(clust=="EV_bedadeti_specific" |
             clust=="EV_mazia_specific" |
             clust=="MA_specific"|
             clust=="MB_specific" ) %>%
    # group_by(clust) %>%
    # summarise(count = n())

    # Join genes who have GO-terms 
    
    left_join(
      GO_terms_global_genes 
    ) %>%
    mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
    filter(!str_detect(GO_terms,"GO")) %>%
    
    # Join genes that showed blastp-based homology
    left_join(
      EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits
    ) %>%
    mutate(
      description = str_replace_na(description, "NA")
    ) %>%
    filter(!str_detect(description,"GO")) %>%
    rename(blastp_annot = description) %>%
    
    # Join protein pfam annotation
    left_join(
      EV_MAB_protinfer %>%
        rename(gene_id = sequence_name) %>%
        filter(!str_detect(predicted_label,"GO:")) %>%
        filter(confidence > 0.9 ) %>%
        separate(gene_id,into = c("gene_id","isoform"),sep = "\\-") 
    ) %>%
    mutate(predicted_label = str_replace_na(predicted_label,"NA")) %>%
    filter(predicted_label !="NA") 
  
  proteinfer_nogo_noblastp_00 %>%
    filter(clust=="EV_mazia_specific") %>%
    distinct(gene_id, predicted_label,description) %>%
    filter(gene_id != "EVBD.1.039879-RA"|
             gene_id != "EVBD.1.046834-RA"|
             gene_id != "EVBD.1.050883-RA" |
             gene_id != "EVMZ.1.003931-RA" |
             gene_id != "EVMZ.1.004895-RA" |
             gene_id != "EVMZ.1.005345-RA" |
             gene_id != "EVMZ.1.014382-RA" |
             gene_id != "EVMZ.1.016368-RA" |
             gene_id != "EVMZ.1.033454-RA" |
             gene_id != "EVMZ.1.035142-RA" |
             gene_id != "Macma4_07_g23090.1" |
             gene_id != "Macma4_10_g13670.1" |
             gene_id != "Mba00_g07130.1"| 
             gene_id != "EVMZ.1.014193" ,
             gene_id != "EVMZ.1.014369", 
             gene_id != "EVMZ.1.014411",
           gene_id != "EVMZ.1.014369",
           gene_id != "EVMZ.1.005015",
           gene_id != "EVMZ.1.015881",
           gene_id != "EVMZ.1.014423",
           gene_id != "EVMZ.1.050165",
           gene_id != "EVMZ.1.049547",
           gene_id != "EVMZ.1.045665",
           gene_id != "EVMZ.1.044045",
           gene_id != "EVMZ.1.041855",
           gene_id != "EVMZ.1.041855",
           gene_id != "EVMZ.1.050103",
           gene_id != "EVMZ.1.030525",
           gene_id != "EVMZ.1.032942",
           gene_id != "EVMZ.1.029637",
           gene_id != "EVMZ.1.014193",
           gene_id != "EVMZ.1.016374",
           gene_id != "EVMZ.1.021573",
           gene_id != "EVMZ.1.004603"
           ) 
  
  
# EV specific protein that lack GO-terms annotation and protein sequences based homology but showed siginficant hits with proteInfer prediction
  
  proteinfer_nogo_noblastp_00 %>%
    filter(clust=="EV_mazia_specific") %>%
    distinct(gene_id,description,predicted_label,confidence) %>%
        filter(gene_id == "EVMZ.1.004895" |
             gene_id == "EVMZ.1.009022" |
             gene_id == "EVMZ.1.014382" |
             gene_id == "EVMZ.1.029956" |
             gene_id == "EVMZ.1.035454"
             ) %>%
    filter(!str_detect(description,"one"),
             !str_detect(description,"Transferases"),
           !str_detect(description,"ester")) %>%
    left_join(
      EV_MAB_prot_fasta %>% 
        separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
    ) %>%
    mutate(isoform = str_replace_na(isoform,"NA")) %>% 
    # filter(isoform !="NA") %>% 
    mutate(seq.name = case_when(!str_detect(isoform,"NA") ~  str_c(gene_id, isoform, sep = "-"),
                                TRUE ~ gene_id)) %>%
    select(seq.name,predicted_label,description,confidence,seq.text) %>%
    write.table(paste0(DATA_DIR, "EV_mazia_nogo_noblastpHomology_butProteinfer.prots.txt"),
                col.names = T,
                row.names = F,
                quote = F,
                sep = '\t')
  

  
  proteinfer_nogo_noblastp_00 %>%
    filter(clust=="EV_mazia_specific") %>%
    distinct(gene_id,description,predicted_label,confidence) %>%
    filter(gene_id == "EVMZ.1.004895" |
             gene_id == "EVMZ.1.009022" |
             gene_id == "EVMZ.1.014382" |
             gene_id == "EVMZ.1.029956" |
             gene_id == "EVMZ.1.035454"
    ) %>%
    filter(!str_detect(description,"one"),
           !str_detect(description,"Transferases"),
           !str_detect(description,"ester")) %>%
    left_join(
      EV_MAB_prot_fasta %>% 
        separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
    ) %>%
    mutate(isoform = str_replace_na(isoform,"NA")) %>% 
    # filter(isoform !="NA") %>% 
    mutate(seq.name = case_when(!str_detect(isoform,"NA") ~  str_c(gene_id, isoform, sep = "-"),
                                TRUE ~ gene_id)) %>%
    select(seq.name,seq.text) %>%
    dat2fasta(paste0(DATA_DIR, "EV_mazia_nogo_noblastpHomology_butProteinfer.prots.fa"))
                
    
    
    
    proteinfer_nogo_noblastp_00 %>%
      # filter(clust=="EV_mazia_specific") %>%
      distinct(gene_id) %>%
      write.table(paste0(DATA_DIR,"EV_MAB_nogo_noblastpHomology_butProteinfer.txt"),
                  col.names = F,
                  row.names = F,
                  quote = F,
                  sep = '\t')
    
    

    proteinfer_nogo_noblastp_00 %>% 
      # filter(clust=="EV_mazia_specific") %>%
      distinct(gene_id) %>%
      left_join(
        EV_MAB_prot_fasta %>% 
          separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
      ) %>%
      mutate(isoform = str_replace_na(isoform,"NA")) %>% 
      # filter(isoform !="NA") %>% 
      mutate(seq.name = case_when(!str_detect(isoform,"NA") ~  str_c(gene_id, isoform, sep = "-"),
                                  TRUE ~ gene_id)) %>%
      # filter(str_detect(isoform,"NA")) %>%
      distinct(seq.name,seq.text) %>%
      filter(str_detect(seq.name,"EVMZ.1.009022"))
      dat2fasta(paste0(DATA_DIR,"EV_MAB_nogo_noblastpHomology_butProteinfer.fa"))
    
    
    proteinfer_nogo_noblastp_00 %>% 
      # filter(clust=="EV_mazia_specific") %>%
      distinct(gene_id) %>%
      left_join(
        EV_MAB_prot_fasta %>% 
          separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
      ) %>%
      mutate(isoform = str_replace_na(isoform,"NA")) %>% 
      # filter(isoform !="NA") %>% 
      mutate(seq.name = case_when(!str_detect(isoform,"NA") ~  str_c(gene_id, isoform, sep = "-"),
                                  TRUE ~ gene_id)) %>%
      filter(str_detect(isoform,"NA")) %>%
      distinct(seq.name,seq.text) %>% 
      dat2fasta(paste0(DATA_DIR,"MAB_nogo_noblastpHomology_butProteinfer.fa"))

    
    proteinfer_nogo_noblastp_00 %>% 
      # filter(clust=="EV_mazia_specific") %>%
      distinct(gene_id,description,predicted_label,confidence) %>%
      left_join(
        EV_MAB_prot_fasta %>% 
          separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
      ) %>%
      mutate(isoform = str_replace_na(isoform,"NA")) %>% 
      # filter(isoform !="NA") %>% 
      mutate(seq.name = case_when(!str_detect(isoform,"NA") ~  str_c(gene_id, isoform, sep = "-"),
                                  TRUE ~ gene_id)) %>%
      filter(str_detect(isoform,"NA")) %>%
      distinct() %>% 
      write.table(paste0(DATA_DIR,"MAB_nogo_noblastpHomology_butProteinfer.txt"),
                  col.names = T,
                  row.names = F,
                  quote = F,
                  sep = "\t")

    
    proteinfer_nogo_noblastp_00 %>%
      filter(str_detect(gene_id,"Macma4_07_g23090"))
  