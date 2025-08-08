
# 
# 
# # Core MAB genes 
# # EV bedadeti gene sets
# EV_bedadeti_specific
# EV_bedadeti_MAB_core
# EV_bedadeti_MA_specific
# EV_bedadeti_MB_specific
# EV_bedadeti_noEVReadSupport_Butsyntenic2MAB
# 
# # EV mazia gene sets
# EV_mazia_specific
# EV_mazia_MAB_core
# EV_mazia_MA_specific
# EV_mazia_MB_specific
# EV_mazia_noEVReadSupport_Butsyntenic2MAB
# 
# # Musa acuminata gene sets
# MA_specific %>% distinct (gene_id) %>% nrow()
# MA_EV_core
# MA_noMAReadSupport_Butsyntenic2EV
# 
# # Musa balbisiana gene sets
# MB_specific
# MB_EV_core
# MB_noMBReadSupport_Butsyntenic2EV


library(phylotools)

EV_MAB_prot_fasta <- c()

for(i in list.files(path = DATA_DIR,pattern = "*.fasta")){
  
  BASENAME <- i %>% 
    str_remove(".fasta") %>%
    str_remove(".prot") 
  
  
  EV_MAB_prot_fasta <- 
    read.fasta(paste0(DATA_DIR,i)) %>%
    mutate(genome = BASENAME) %>%
    rbind(EV_MAB_prot_fasta)
}


EV_MAB_prot_fasta %>%
  # select(seq.name) %>%
  tail(n=1)

# Bedadeti specific 

EV_bedadeti_specific %>%
  distinct() %>%
  left_join(
    # Add full gene id
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "bedadeti") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id)) %>%
  select(seq.name, gene_id) %>% 
  distinct() %>%
  left_join(
    # join with protein sequence 
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_bedadeti") 
  ) %>%
  mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  filter(seq.text != "NA") %>%
  distinct() %>%
  select(-gene_id,-genome) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_bedadeti_specific.01.fa"))

# Mazia specific 

EV_mazia_specific %>%
  distinct() %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "mazia") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id)) %>%
  select(seq.name, gene_id) %>% 
  distinct() %>% 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_mazia") 
  ) %>%
  mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  filter(seq.text != "NA") %>%
  distinct() %>%
  select(-gene_id,-genome) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_mazia_specific.01.fa"))


# MA specific 
MA_specific %>%
  distinct() %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_acuminata") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id)) %>%
  select(seq.name, gene_id) %>% 
  distinct() %>%  
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ac") 
  ) %>%
  mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  filter(seq.text != "NA") %>%
  distinct() %>%
  select(-gene_id,-genome) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MA_specific.01.fa"))


# MB specifc 

MB_specific %>%
  distinct() %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_balbisiana") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id)) %>%
  select(seq.name, gene_id) %>% 
  distinct() %>%  
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ba") 
  ) %>%
  mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  filter(seq.text != "NA") %>%
  distinct() %>%
  select(-gene_id,-genome) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_specific.01.fa"))


# Species/genotype specific genes who do not have go terms annotation
# bedadeti specific proteins 
## get the orginal gene id 
EV_bedadeti_specific %>%
  filter(!str_detect(GO_terms,"GO")) %>% 
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "bedadeti") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_bedadeti") %>% 
      select(-genome)
  ) %>% 
  select(seq.name,seq.text) %>%
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_bedadeti_specific_nogo.syn.01.fa"))


#  Mazia specific proteins 
## get the orginal gene id 
EV_mazia_specific %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "mazia") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_mazia") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>% 
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_mazia_specific_nogo.syn.01.fa"))


# for MA specific proteins 
## get the orginal gene id 
MA_specific %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_acuminata") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ac") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","musa_ac_specific_nogo.syn.01.fa"))


# for MB proteins 
## get the orginal gene id 
MB_specific %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_balbisiana") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ba") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct(seq.name) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","musa_ba_specific_nogo.01.fa"))


### genes that lack self-read support (<0.25 coverage)

EV_bedadeti_noEVReadSupport_Butsyntenic2MAB %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "bedadeti") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_bedadeti") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct() %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_bedadeti_noEVReadSupport_Butsyntenic2MAB_nogo.01.fa"))

EV_MAB_prot_fasta %>%
  filter(str_detect(seq.name,"EVMZ.1.037942"))

EV_mazia_noEVReadSupport_Butsyntenic2MAB %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "mazia") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>% 
  select(seq.name) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "EV_mazia") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_mazia_noEVReadSupport_Butsyntenic2MAB_nogo.01.fa"))


MA_noMAReadSupport_Butsyntenic2EV %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_acuminata") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ac") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MA_noMAReadSupport_Butsyntenic2EV_nogo.01.fa"))



MB_noMBReadSupport_Butsyntenic2EV %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  left_join(
    GO_terms_EV_AED25_MAB_desc.v1 %>%
      filter(genome == "Musa_balbisiana") %>%
      mutate(gene_id=str_extract(seq.name,"\\w+.\\d+.\\d+"),
             GO_terms = str_replace_na(GO_terms,"NA")) %>%
      select(seq.name,gene_id,GO_terms)) %>%
  ## bind protein sequences 
  left_join(
    EV_MAB_prot_fasta %>%
      filter(genome == "musa_ba") %>% 
      select(-genome)
  ) %>%
  select(seq.name,seq.text) %>%
  distinct() %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_noMBReadSupport_Butsyntenic2EV_nogo.01.fa"))


# Get CDS of EV and MAB specific genes
EV_MAB_reads_coverage_geneID %>%
  filter(ref_genome=="musa_acuminata") %>%
  head()


#Species / genotype specific genes 

EV_bedadeti_specific %>%
  rbind(EV_bedadeti_noEVReadSupport_Butsyntenic2MAB) %>%
  select(-GO_terms) %>%
  distinct() %>%
  mutate(cat = "EV_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      rename(gene_id = Name)
  ) %>%
  select(seqid,start,end,gene_id) %>%
  select(seqid,start,end) %>%
  distinct() %>%
  write.table("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir/EV_bedadeti_spec_plus_noreadSup.01.bed",
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")

EV_mazia_specific %>%
  rbind(EV_mazia_noEVReadSupport_Butsyntenic2MAB) %>%
  select(-GO_terms) %>%
  distinct() %>%
  mutate(cat = "EV_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      rename(gene_id = Name)
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(seqid,start,end) %>%
  distinct() %>%
  write.table("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir/EV_mazia_spec_plus_noreadSup.01.bed",
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


MA_specific %>%
  rbind(MA_noMAReadSupport_Butsyntenic2EV) %>%
  select(-GO_terms) %>%
  distinct() %>%
  mutate(cat = "MA_specific") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      filter(ref_genome=="musa_acuminata") %>%
      rename(gene_id = Name) %>% 
      mutate(gene_id = str_c(gene_id,".1"))
  ) %>% 
  # distinct(gene_id) %>%
  # nrow()
  select(seqid,start,end) %>%
  distinct() %>% 
  write.table("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir/MA_spec_plus_noreadSup.01.bed",
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")

MB_specific %>%
  rbind(MB_noMBReadSupport_Butsyntenic2EV) %>%
  select(-GO_terms) %>% 
  distinct() %>%
  mutate(cat = "MB_specific") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      filter(ref_genome=="musa_balbisiana") %>%
      rename(gene_id = Name) %>% 
      mutate(gene_id = str_c(gene_id,".1"))
  ) %>% 
  # distinct(gene_id) %>%
  # nrow()
  select(seqid,start,end) %>%
  distinct() %>% 
  write.table("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir/MB_spec_plus_noreadSup.01.bed",
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")

