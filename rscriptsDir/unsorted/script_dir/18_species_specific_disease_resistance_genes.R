#EV-MAB genes 

EV_MAB_prot_fasta

EV_MAB_reads_coverage_geneID %>%
  filter(str_detect(Name,"EVBD.1.033364")|
           str_detect(Name,"EVBD.1.051740")|
           str_detect(Name,"EVMZ.1.034720")) %>%
  head()

EV_MAB_reads_coverage_geneID_synteny %>%
  head()

EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  filter(str_detect(gene_id,"EVBD.1.033364")|
           str_detect(gene_id,"EVBD.1.051740")|
           str_detect(gene_id,"EVMZ.1.034720")) %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>%
  write.table(paste0(DATA_DIR,"EV_specific_Rgenes_nogo_lack_homology.bed"),
              col.names = F,
              row.names = F,
              quote= F,
              sep = "\t")


EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  filter(str_detect(gene_id,"EVBD.1.033364")|
           str_detect(gene_id,"EVBD.1.051740")|
           str_detect(gene_id,"EVMZ.1.034720")) %>%
  select(seqid,start,end,gene_id) %>%
  left_join(
    EV_MAB_prot_fasta %>% 
      separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
  ) %>%
  mutate(isoform = str_replace_na(isoform,"NA"),
         gene_id = str_c(gene_id,isoform, sep="-")) %>% 
  select(-genome) %>%
  distinct()  %>%
  write.table(paste0(DATA_DIR,"EV_specific_Rgenes_nogo_lack_homology.txt"),
              col.names = F,
              row.names = F,
              quote= F,
              sep = "\t")

EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  filter(str_detect(gene_id,"EVBD.1.033364")|
           str_detect(gene_id,"EVBD.1.051740")|
           str_detect(gene_id,"EVMZ.1.034720")) %>%
  select(seqid,start,end,gene_id) %>%
  left_join(
    EV_MAB_prot_fasta %>% 
      separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
  ) %>%
  mutate(isoform = str_replace_na(isoform,"NA"),
         gene_id = str_c(gene_id,isoform, sep="-")) %>% 
  rename(seq.name = gene_id) %>%
  distinct(seq.name,seq.text)  %>%
  dat2fasta(paste0(DATA_DIR,"EV_specific_Rgenes_nogo_lack_homology.fa"))




EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  mutate(
    description = str_replace(description, "hypersensitive induced reaction","Hypersensitive-induced response protein"),
    description = str_replace(description, "disease resistance protein","Disease resistance protein")) %>%
  
  filter(str_detect(description,"Zinc finger proteins")) %>%
  distinct(gene_id) %>%
  left_join(
    EV_MAB_prot_fasta %>% 
      separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
  ) %>%
  mutate(isoform = str_replace_na(isoform,"NA"),
         gene_id = str_c(gene_id,isoform, sep="-")) %>% 
  rename(seq.name = gene_id) %>%
  dat2fasta(paste0(DATA_DIR,"EV_specific_nogo_blasthomology_finger_proteins.fa"))






EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  # extract genes that showed blastp similary finger proteins
  filter(str_detect(gene_id,"EVMZ.1.007075")|
           str_detect(gene_id,"EVMZ.1.022832")|
           str_detect(gene_id,"Macma4_03_g23940.1")|
           str_detect(gene_id,"Macma4_04_g03220.1")|
           str_detect(gene_id,"Macma4_04_g18850.1")|
           str_detect(gene_id,"Macma4_11_g00430.1")) %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>%
  # Join with protein sequences
  left_join(
    EV_MAB_prot_fasta %>% 
      separate(seq.name,into = c("gene_id","isoform"),sep = "\\-") 
  ) %>%
  mutate(isoform = str_replace_na(isoform,"NA"),
         gene_id = str_c(gene_id,isoform, sep="-")) %>%
  select(-isoform) %>%
  write.table(paste0(DATA_DIR,"EV_specific_nogo_blasthomology_finger_proteins.txt"),
              col.names = T,
              row.names = F,
              quote= F,
              sep = "\t")






