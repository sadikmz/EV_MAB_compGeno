# 

read.fasta(paste0("Z:/sadik/evmab_dir/","EV_MAB_novel_genes_dmdb/EV_MAB.novel_genes.v2.fa")) %>%
rbind(
read.fasta(paste0("Z:/sadik/evmab_dir/","phylogene_dir/zingiberaceae_chpt_mito.pro.fa"))) %>%
  rownames_to_column() %>%
  mutate(seqname = case_when(str_detect(seq.name,"^EVMZ") ~ "EVMZ",
                             str_detect(seq.name,"^EVBD") ~ "EVBD",
                             str_detect(seq.name,"^Ma") ~ "Ma",
                             str_detect(seq.name,"^Mb") ~ "Mb" ,
                             TRUE ~ "rt")) %>%
  mutate(seq_name = str_c(rowname,seqname,sep = "-")) %>%
  select(seq_name,seq.text) %>%
  rename(seq.name = seq_name) %>%nrow()
  dat2fasta(paste0("Z:/sadik/evmab_dir/","phylogene_dir/EV_MAB.novel_vs_zingibrales_plastome.fa"))
  
  
# Generate phylogenetic tree for genes whose alphafold proteins structure confidence > 50 and < 50s