EV_MAB_reads_coverage_geneID_GO_terms_synteny_uce <- 
  
  # All core EV-MAB genes 
  
  EV_bedadeti_MAB_core %>%
  rbind(
    EV_mazia_MAB_core,
    MA_EV_core,
    MB_EV_core) %>%
  mutate(cat = "core") %>%
  
  # Extract and subset with genes with 100% lastz alignment and have >=70 % interspecies reads coverage. 
  # These core genes with 100% lastz alignment +70 interspecies reads coverage were considered as UCE (ultra conserved elements)
  left_join(
EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%

  mutate(reads_genome = str_replace(reads_genome,"EV genotypes", "EV"),
         reads_genome = str_replace(reads_genome,"A-sub genome banana genotypes", "MA"),
         reads_genome = str_replace(reads_genome,"B-sub genome banana genotypes", "MB")
  ) %>% 
  # head()
  # distinct(reads_used)
  dplyr::filter(reads_used != "qfiltered") %>%
  dplyr::filter(gene_repeat == "gene") %>%
  mutate(syn_cov = str_replace_na(syn_cov,"NA"),
         synteny = str_replace_na(synteny,"NA")) %>%
  dplyr::select(gene_id, ref_genome, query, reads_used, reads_genome,synteny,GO_terms,fraction_overlap,syn_cov) %>%
  filter(fraction_overlap  >=0.70) %>%
  filter(syn_cov == 1) %>%
  select(-synteny ) %>%
  as.data.frame () %>% 
  distinct() ) %>%
  mutate(cat = str_replace_na(cat, "NA"),
         ref_genome = str_replace_na(ref_genome,"NA")) %>%
  filter(cat != "NA") %>%
  filter(ref_genome !="NA") 



EV_MAB_reads_coverage_geneID_GO_terms_synteny_uce %>% 
  select(ref_genome,gene_id) %>% 
  distinct() %>%
    group_by(ref_genome) %>%
    summarise(count = n())

EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  
  mutate(reads_genome = str_replace(reads_genome,"EV genotypes", "EV"),
         reads_genome = str_replace(reads_genome,"A-sub genome banana genotypes", "MA"),
         reads_genome = str_replace(reads_genome,"B-sub genome banana genotypes", "MB")
  ) %>% 
  dplyr::filter(reads_used != "qfiltered") %>%
  dplyr::filter(gene_repeat == "gene") %>%
  filter(GO_terms == "GO:0015074") %>%
  filter(synteny!="self_reads") %>%
  filter(synteny=="non_syntenic") %>%
  select(genome,gene_id,fraction_overlap ) %>%
  distinct() 
  summarise(count = n())

  
  library(tidyverse)
  library(phylotools)
  
  DATA_DIR="/home/lifesci/lfrwtp/data/genomes/"


  for(i in list.files(path=DATA_DIR, pattern="spec_plus_noreadSup.fna")){
    
    BASENAME=basename(i) %>% str_remove("_spec_plus_noreadSup.fna")
    
    read.fasta(paste0(DATA_DIR,i)) %>%
      rownames_to_column() %>%
      mutate(seqname = str_c(BASENAME,"_",rowname)) %>%
      select(seqname,seq.text) %>%
      rename(seq.name = seqname) %>%
      dat2fasta(paste0(DATA_DIR, BASENAME,"_spec_plus_noreadSup.renamed.fna"))
    
    read.fasta(paste0(DATA_DIR,i)) %>%
      rownames_to_column() %>%
      mutate(seqname = str_c(BASENAME,"_",rowname)) %>%
      select(seqname,seq.name) %>%
      write.table(paste0(DATA_DIR, BASENAME,"_spec_plus_noreadSup.renamed.txt"),
                  col.names = T, row.names = F, quote = F, sep = '\t')
  }

