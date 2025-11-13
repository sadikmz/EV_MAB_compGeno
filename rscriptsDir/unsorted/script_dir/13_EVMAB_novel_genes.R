

# EVMAB that lack go-terms annotation as well as blastp against ncbi nr, uniprot and at and os
EV_MAB_specific_genes_nogo_noblastp %>%
  head()
# EV_MBA gene clusters 
EV_MAB_custered_genes.v1 %>%
  head()

# Interproscan 
EV_MAB_interprocan %>%
  head()

#Proteinfer 
EV_MAB_protinfer

# eggnog 

eggnog_MAB <-
  read.delim(paste0(DATA_DIR,"eggnog_mapper_musa_ac.out/eggNOG.emapper.annotations"),skip = 4, sep = '\t') %>%
  rbind(
    read.delim(paste0(DATA_DIR,"eggnog_mapper_musa_ac.out/eggNOG.emapper.annotations"),skip = 4, sep = '\t'),
    read.delim(paste0(DATA_DIR,"eggnog_mapper_bedadeti_AED25.out/eggNOG.emapper.annotations"),skip = 4, sep = '\t'),
    read.delim(paste0(DATA_DIR,"eggnog_mapper_mazia_AED25.out/eggNOG.emapper.annotations"),skip = 4, sep = '\t')
  ) %>%
  rename(gene_id = X.query) %>%
  mutate(count_sep = str_count(GOs,",")) %>%
  separate(gene_id,c("gene_id","isomer"), sep = "-") %>%
  distinct

eggnog_MAB %>%
  head()


# Extract EV_MAB specific genes that lack go temrs annnotation and do not have hits against ncbi nr, uniprot, at_os proteins as well as eggnog annotation, interproscan and proteinfer search



# EV_MBA gene clusters 
EV_MAB_custered_genes.v2 <-
# EV_MAB_custered_genes.v1 %>%
# join with genes lack go temrs annnotation and do not have hits against ncbi nr, uniprot, at_os proteins

EV_MAB_specific_genes_nogo_noblastp %>%
  mutate (annot_no_blastp = "no_blastp") %>%
  rename(cds_start = start,
         cds_end = end) %>%
  left_join(
    EV_MAB_custered_genes.v1) %>%
  # Join interproscan result of EV_MAB proteom
  
    left_join(EV_MAB_interprocan %>%
              rename(gene_id=seq.name) %>%
                mutate(annot_intpro = "interprocan") %>%
  rename(seq_mds = V2,
         seq_len = V3,
         analysis = V4,
         signature_accession = V5,
         signature_description = V6,
         intpro_start = V7,
         intpro_end = V8,
         score = V9,
         status = V10,
         date = V11,
         interproscan_annotaton = V12,
         interproscan_description = V13,
         GO_annotation=V14
  ) %>%
  
  select(-seq_mds,-date) %>%
  mutate(analysis = str_replace_na(analysis,"NA"),
         signature_accession = str_replace_na(signature_accession,"NA"),
         signature_description = str_replace_na(signature_description,"NA"),
         interproscan_annotaton = str_replace_na(interproscan_annotaton,"NA"),
         interproscan_description = str_replace_na(interproscan_description,"NA")) %>%
  # filter(analysis != "NA") %>%
  filter(interproscan_description == "-") %>%
  filter(signature_description !="consensus disorder prediction",
         signature_description != "Coil",
         signature_description != "-")) %>%
  
  
  # Join Proteinfer result of EV_MAB proteom
  
  left_join(
    EV_MAB_protinfer %>%
      filter(!str_detect(predicted_label,"GO:")) %>%
      filter(confidence > 0.9 ) %>%
      separate(sequence_name,into=c("gene_id","isomer"),sep="-") %>%
      mutate(predicted_label  = str_replace_na(predicted_label, "NA")) %>%
      filter(predicted_label !="NA") %>%
      rename(protinfer_description = description) %>%
      select(-isomer) %>%
      mutate(annot_protenfer = "proinfer")
  ) %>%

  # Join eggnog result of EV_MAB proteom
  
  left_join(
    eggnog_MAB %>%
      select(gene_id,eggNOG_OGs,Description, Preferred_name) %>%
      mutate(annot_eggnog = "eggnog")
    
  ) %>%
  mutate(annot_eggnog = str_replace_na(annot_eggnog,"NA"),
         annot_protenfer = str_replace_na(annot_protenfer,"NA"),
         annot_intpro = str_replace_na(annot_intpro,"NA"),
         mapping_reads_source = str_replace_na(mapping_reads_source,"NA"),
         predicted_label = str_replace_na(predicted_label,"NA"),) 



# 
EV_MAB_specific_genes_nogo_noblastp %>%
  filter(mapping_reads_source == "cross_species_reads") %>%
  filter( frac_ovrerlap_syn < 0.245) %>%
  filter(genome== "mazia") %>%
distinct(gene_id) %>% nrow()


# EV_MAB nogo_noblastp genes that showed proinfer, interproscan and eggnog annotation
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


# Final genes novel EV_MAB genes list
EV_MAB_custered_genes.v3 <- 
EV_MAB_custered_genes.v2 %>% 
  filter(mapping_reads_source == "cross_species_reads") %>%
  # filter( query != "MS_LongReads",
  #           query != "MT_LongReads") %>%
  filter( frac_ovrerlap_syn < 0.245) %>%
  select(seqid, cds_start, cds_end, genome, gene_id, predicted_label, protinfer_description,eggNOG_OGs,Preferred_name) %>%
  full_join(
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
  filter(protinfer_description != "NA") %>%
  distinct(gene_id) %>%
  mutate(annot = "proteinfer_annot")) %>%
  mutate(annot = str_replace_na(annot, "NA")) %>%
  filter(annot != "proteinfer_annot") %>%
    distinct(seqid,cds_start,cds_end,gene_id,genome)



EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "mazia") %>%
  select(-gene_id,-genome) %>%
  distinct() %>%
    write.table(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.bed"),
            col.names = F, 
            row.names = F,
            quote = F,
            sep = "\t")

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "bedadeti") %>%
  select(-gene_id,-genome) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_acuminata") %>%
  select(-gene_id,-genome) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.genes.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  select(-gene_id,-genome) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.genes.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


# EV mazia genes that lack read support from EV genotypes but show synteny with MAB

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "mazia") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  select(seq.name,seq.text) %>%
dat2fasta(paste0("Z:/sadik/evmab_dir/","mazia_specific.novel.genes.fa"))


EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "bedadeti") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  select(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","bedadeti_specific.novel.genes.fa"))


EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_acuminata") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  select(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MA_specific.novel.genes.fa"))

EV_MAB_custered_genes.v3 %>%
  # distinct(genome)
  filter(genome == "Musa_balbisiana") %>%
  distinct(gene_id) %>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  select(seq.name,seq.text) %>% 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_specific.novel.genes.fa"))



