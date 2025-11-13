# com bined genes list of EV and MAB gene lists
EV_MAB_custered_genes.v4 %>%
  distinct(clust)

# Extract CDS coords and save to disk 
EV_MAB_custered_genes.v4 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      select(seqid,start,end,gene_id,ref_genome) %>%
  distinct()) %>% 
  filter(ref_genome == "ev_bedadeti") %>%
  filter(clust == "EV_bedadeti_specific") %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>% head()
  write.table(paste0(venndiagOut,"EV_bedadeti_specific_CDS.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")

# Mazia
EV_MAB_custered_genes.v4 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      select(seqid,start,end,gene_id,ref_genome) %>%
      distinct()) %>%
  filter(ref_genome == "ev_mazia") %>%
  filter(clust == "EV_mazia_specific") %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>%
  # distinct(gene_id) %>% nrow()
# head() 
  write.table(paste0(venndiagOut,"EV_mazia_specific_CDS.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")

# MA
EV_MAB_custered_genes.v4 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      select(seqid,start,end,gene_id,ref_genome) %>%
      distinct()) %>%
  # distinct(ref_genome)
  filter(ref_genome == "musa_acuminata") %>%
  filter(clust == "MA_specific") %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>% 
  # nrow()
  # head()
  write.table(paste0(venndiagOut,"MA_specific_CDS.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")

# MB

EV_MAB_custered_genes.v4 %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      select(seqid,start,end,gene_id,ref_genome) %>%
      distinct()) %>%
  # distinct(ref_genome)
  filter(ref_genome == "musa_balbisiana") %>%
  filter(clust == "MB_specific") %>%
  select(seqid,start,end,gene_id) %>%
  distinct() %>%
  # nrow()
  # head()
  write.table(paste0(venndiagOut,"MB_specific_CDS.bed"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")

# rename fasta seq

library(phylotools)

fastaList=list.files(path = "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/",
             pattern = "CDS.fna")

# fastaIn=paste0("Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/MA_specific_CDS.fna")

for (FASTA in fastaList){
  
  # Basename 
  
  BASENAME=basename(FASTA) %>% str_remove(".fna")
  
  # read and rename sequence names 
  read.fasta(paste0("Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/", FASTA)) %>%
    rownames_to_column() %>%
    # select(-seq.text) %>%
    mutate(seq_name = str_c("seq_",rowname),
           seq_len = nchar(seq.text)) %>%
    select(seq_name,seq.name,seq_len) %>%
    write.table(paste0("Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/", BASENAME, ".00.txt"),
                row.names = F,
                col.names = T,
                quote = F,
                sep = "\t")
  
  
  
  read.fasta(paste0("Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/", FASTA)) %>%
    rownames_to_column() %>%
    # select(-seq.text) %>%
    mutate(seq_name = str_c("seq_",rowname)) %>%
    select(seq_name,seq.text) %>%
    rename(seq.name = seq_name) %>%
    dat2fasta(paste0("Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/", BASENAME, ".00.fna"))
  
}


# After interspecies lastz alignment and coverage alignment added as lastz_cov

EV_MAB_custered_genes.v5 <-
  # MA
  read.delim(
    paste0(
      "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/",
      "MA_specific_CDS.00.txt"
    ),
    header = T
  ) %>%
  # head() %>%
  left_join(
    read.delim(
      paste0(
        "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/mamb/",
        "ma_lastz.bed"
      ),
      header = F
    ) %>%
      mutate(cov = round(V4 / (V3 - V2 + 1), 2)) %>%
      rename(seq_name = V1, lastz_cov = cov)
  ) %>%
  separate(seq.name, into = c("seqid", "coord"), sep = ":") %>%
  separate(coord, into = c("start", "end"), sep = "-") %>%
  select(-V2, -V3, -V4, -seq_len, -seq_name) %>%
  distinct() %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  left_join(EV_MAB_reads_coverage_geneID_GO_terms) %>%
  select(gene_id, lastz_cov) %>%
  distinct() %>%
  mutate(clust = "MA_specific") %>% 
  
  # bind
  rbind(
    # MB
    read.delim(
      paste0(
        "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/",
        "MB_specific_CDS.00.txt"
      ),
      header = T
    ) %>%
      left_join(
        read.delim(
          paste0(
            "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/mbma/",
            "mb_lastz.bed"
          ),
          header = F
        ) %>%
          mutate(cov = round(V4 / (V3 - V2 + 1), 2)) %>%
          rename(seq_name = V1, lastz_cov = cov)
      ) %>%
      separate(seq.name, into = c("seqid", "coord"), sep = ":") %>%
      separate(coord, into = c("start", "end"), sep = "-") %>%
      select(-V2, -V3, -V4, -seq_len, -seq_name) %>%
      distinct() %>%
      mutate(start = as.numeric(start), end = as.numeric(end)) %>%
      left_join(EV_MAB_reads_coverage_geneID_GO_terms) %>%
      select(gene_id, lastz_cov) %>%
      mutate(clust = "MB_specific") %>%
      distinct(),
    
    
    
    # EV mazia
    read.delim(
      paste0(
        "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/",
        "EV_mazia_specific_CDS.00.txt"
      ),
      header = T
    ) %>%
      left_join(
        read.delim(
          paste0(
            "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/evmz/",
            "mazia_lastz.bed"
          ),
          header = F
        ) %>%
          mutate(cov = round(V4 / (V3 - V2 + 1), 2)) %>%
          rename(seq_name = V1, lastz_cov = cov)
      ) %>%
      separate(seq.name, into = c("seqid", "coord"), sep = ":") %>%
      separate(coord, into = c("start", "end"), sep = "-") %>%
      select(-V2, -V3, -V4, -seq_len, -seq_name) %>%
      # distinct() %>%
      mutate(start = as.numeric(start), end = as.numeric(end)) %>%
      left_join(EV_MAB_reads_coverage_geneID_GO_terms) %>%
      select(gene_id, lastz_cov) %>%
      mutate(clust = "EV_mazia_specific") %>%
      distinct(),
    
    
    # EV bedadeti
    read.delim(
      paste0(
        "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/fasta_dir/",
        "EV_bedadeti_specific_CDS.00.txt"
      ),
      header = T
    ) %>%
      left_join(
        read.delim(
          paste0(
            "Z:/sadik/evmab_dir/data_dir/synteny_dir/lastz_species_specific_dir/evbd/",
            "bedadeti_lastz.bed"
          ),
          header = F
        ) %>%
          mutate(cov = round(V4 / (V3 - V2 + 1), 2)) %>%
          rename(seq_name = V1, lastz_cov = cov)
      ) %>%
      separate(seq.name, into = c("seqid", "coord"), sep = ":") %>%
      separate(coord, into = c("start", "end"), sep = "-") %>%
      select(-V2, -V3, -V4, -seq_len, -seq_name) %>%
      distinct() %>%
      mutate(start = as.numeric(start), end = as.numeric(end)) %>%
      left_join(EV_MAB_reads_coverage_geneID_GO_terms) %>%
      select(gene_id, lastz_cov) %>%
      mutate(clust = "EV_bedadeti_specific") %>%
      distinct()
  ) %>%
  distinct()


EV_MAB_custered_genes.v5 %>%
  # group_by(clust) %>%
  # summarise(count = n())
  head()

# EV_MAB_custered_genes.v4 is an updated version but interspecies lastz or interaspeces lastz alignment was done on the prevouse output lastz alignmen containing
# 3 MA, 27 EV mazia and 20 EV mazia genes that showed 24.4%<cov<25% lastz alignment coverage. 
# Here these marginaly aligning genes were filtered out for species specific list. 

EV_MAB_custered_genes.v5.1 <- 
EV_MAB_custered_genes.v4 %>%
  mutate(source = "marginal_lastz_filtered") %>%
  left_join(EV_MAB_custered_genes.v5) %>% 
  mutate(source = str_replace_na(source, "NA")) %>%
  distinct(gene_id,lastz_cov,clust) 
  
EV_MAB_reads_coverage_geneID_GO_terms %>%
  distinct(gene_id,GO_terms) %>%
  head()


EV_MAB_custered_genes.v5.1 %>%
  # group_by(clust) %>%
  # summarise(count = n())
  head()
