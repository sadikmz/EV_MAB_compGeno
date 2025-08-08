# load  from 14_EVMAB_novel_genes_blastp.R

EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
  distinct(qseqid) %>%
  arrange(qseqid) %>%
  filter(str_detect(qseqid,"EV")) %>%
  distinct(qseqid) %>%
  separate(qseqid,into=c("qseqid","isomer"), sep="-")


# Mazia
EV_MAB_custered_genes.v3 %>%
  filter(genome == "mazia") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  # select(-gene_id,-genome) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer != "NA") %>%
  select(-gene_id,-genome,-isomer) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.with_blastp_hits.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


EV_MAB_custered_genes.v3 %>%
  filter(genome == "mazia") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  # select(-gene_id,-genome) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer == "NA") %>%
  select(-gene_id,-genome,-isomer) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.refined_no_blastp.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")

# Bedadeti 

EV_MAB_custered_genes.v3 %>%
  filter(genome == "bedadeti") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  # select(-gene_id,-genome) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer != "NA") %>%
  select(-gene_id,-genome,-isomer) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.with_blastp_hits.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


EV_MAB_custered_genes.v3 %>%
  filter(genome == "bedadeti") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  # select(-gene_id,-genome) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer == "NA") %>%
  select(-gene_id,-genome,-isomer) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.refined_no_blastp.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")



## MA

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_acuminata") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type != "NA") %>%
  select(-gene_id,-genome,-type) %>%
  write.table(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.genes_with_blastp_hits.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")



EV_MAB_custered_genes.v3 %>%
  filter(genome == "Musa_acuminata") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type == "NA") %>%
  select(-gene_id,-genome,-type) %>%
  write.table(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.refined_no_blasts.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")



#MB

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type != "NA") %>%
  select(-gene_id,-genome,-type) %>%
  # distinct() 
  write.table(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.genes_with_blastp_hits.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")



EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type == "NA") %>%
  select(-gene_id,-genome,-type) %>%
  # distinct() 
  write.table(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.refined_no_blasts.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


# Refined EV-MAB novel genes 
EV_MAB_custered_genes.v3 %>%
  filter(genome == "mazia") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer != "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.with_blastp_hits.fa"))

EV_MAB_custered_genes.v3 %>%
  filter(genome == "mazia") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer == "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.refined_no_blastp.fa"))

EV_MAB_custered_genes.v3 %>%
  filter(genome == "bedadeti") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer != "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.with_blastp_hits.fa"))

EV_MAB_custered_genes.v3 %>%
  filter(genome == "bedadeti") %>%
  left_join(EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
              distinct(qseqid) %>%
              arrange(qseqid) %>%
              filter(str_detect(qseqid,"EV")) %>%
              distinct(qseqid) %>%
              separate(qseqid,into=c("qseqid","isomer"), sep="-") %>%
              rename(gene_id = qseqid)) %>%
  mutate(isomer = str_replace_na(isomer,"NA")) %>%
  filter(isomer == "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.refined_no_blastp.fa"))

## MA

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_acuminata") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type != "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.genes.with_blastp_hits.fa"))

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_acuminata") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type == "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.genes.refined_no_blastp.fa"))

  ## MB

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type != "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.genes.with_blastp_hits.fa"))

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  left_join(
    
    EV_MAB.novel_genes_ncbi_nr_uniprot_atos_blastp %>%
      distinct(qseqid) %>%
      arrange(qseqid) %>%
      rename(gene_id = qseqid) %>%
      mutate(type ="blastp")
  ) %>% 
  mutate(type = str_replace_na(type,"NA")) %>%
  filter(type == "NA") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  distinct(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.genes.refined_no_blastp.fa"))



