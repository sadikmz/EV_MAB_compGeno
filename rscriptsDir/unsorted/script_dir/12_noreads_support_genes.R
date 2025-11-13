# MB genes that lack read support from MB genotypes but show synteny with EV
MB_noMBReadSupport_Butsyntenic2EV %>%
  distinct(gene_id)%>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta) %>%
  select(seq.name,seq.text) %>%
  dat2fasta(paste0("Z:/sadik/evmab_dir/","MB_noMBReadSupport_Butsyntenic2EV.fa"))



# Extract coding sequence 
MB_noMBReadSupport_Butsyntenic2EV %>%
  distinct(gene_id)%>%
  left_join(EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1) %>%
  filter(reads_used != "qfiltered") %>%
  select(seqid,start,end,gene_id,frac_ovrerlap_syn,syn_cov,GO_terms) %>%
  mutate (frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) %>%
  distinct() %>%
  write.table(paste0("Z:/sadik/evmab_dir/","MB_noMBReadSupport_Butsyntenic2EV.geneID_fracOvlp_syn.bed"),
              col.names = F, 
              row.names = F,
              quote = F,
              sep = "\t")


# EV mazia genes that lack read support from EV genotypes but show synteny with MAB

EV_mazia_noEVReadSupport_Butsyntenic2MAB %>%
  distinct(gene_id)%>%
  rename(seq.name = gene_id) %>%
  left_join(EV_MAB_prot_fasta %>%
              separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
  mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
  select(seq.name,seq.text) 
  dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_mazia_noEVReadSupport_Butsyntenic2MAB.fa"))

  EV_mazia_noEVReadSupport_Butsyntenic2MAB %>%
    distinct(gene_id)%>%
    left_join(EV_MAB_reads_coverage_geneID_GO_terms_synteny.v1) %>%
    filter(reads_used != "qfiltered") %>%
    rename(seq.name = gene_id) %>% 
    select(-genome) %>%
    left_join(EV_MAB_prot_fasta %>%
                separate(seq.name,c("seq.name","isomer"), sep="-")) %>%
    mutate(seq.name = str_c(seq.name,isomer,sep = "-")) %>%
    select(seqid,start,end,seq.name,frac_ovrerlap_syn,syn_cov,GO_terms) %>%
    mutate (frac_ovrerlap_syn = round(frac_ovrerlap_syn,2)) %>%
    distinct() %>%
    dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_mazia_noEVReadSupport_Butsyntenic2MAB.geneID_fracOvrlp_syn.GO.bed"))
  
  
  
    EV_MAB_prot_fasta %>%
    select(seq.name,seq.text) %>%
    dat2fasta(paste0("Z:/sadik/evmab_dir/","EV_MAB.prot..fa"))
  