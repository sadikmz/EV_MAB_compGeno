# extract glycerol-3-phosphate gene 
blast_out.1 %>%
  filter(
    #description=="carbohydrate metabolic process" |
    str_detect(stitle_bp_annot, "ransposable") |
    # str_detect(description, "ranscription")
    str_detect(description, "DNA integration")
    # str_detect(description, "glutathione metabolic process")
  ) %>%
  select(seqid, start, end,seq.name) %>%
  filter(!str_detect(seq.name,"EVMZ")) %>%
  distinct(seqid,start,end) %>%
  write.table(paste0(result_dir,"MA_specific_DNA_intergration_transposable.cds.bed"),
              col.names = F,
              row.name = F,
              quote = F,
              sep = "\t")



blast_out.1 %>%
  filter(
    #description=="carbohydrate metabolic process" |
    str_detect(stitle_bp_annot, "ransposable") |
      # str_detect(description, "ranscription")
      str_detect(description, "DNA integration")
    # str_detect(description, "glutathione metabolic process")
  ) %>%
  select(seqid, start, end,seq.name) %>%
  filter(str_detect(seq.name,"EVMZ")) %>%
  distinct(seqid,start,end) %>%
  write.table(paste0(result_dir,"EVMZ_specific_DNA_intergration_transposable.cds.bed"),
              col.names = F,
              row.name = F,
              quote = F,
              sep = "\t")

# write predicted proteins 
blast_out.1 %>%
  filter(
    #description=="carbohydrate metabolic process" |
    str_detect(stitle_bp_annot, "ransposable") |
      # str_detect(description, "ranscription")
      str_detect(description, "DNA integration")
    # str_detect(description, "glutathione metabolic process")
  ) %>%
  select(seqid, start, end,seq.name) %>%
  filter(!str_detect(seq.name,"EVMZ")) %>%
  distinct() %>%
  left_join(
    EV_MAB_prot_fasta %>%
      filter(!str_detect(seq.name,"EV")) %>%
      select(-genome)
  ) %>%
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>%
  write.table(paste0(result_dir, "MAB_specific_transposable.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)




blast_out.1 %>%
filter(
  #description=="carbohydrate metabolic process" |
  str_detect(stitle_bp_annot, "ransposable") |
    # str_detect(description, "ranscription")
    str_detect(description, "DNA integration")
  # str_detect(description, "glutathione metabolic process")
) %>%
  select(seqid, start, end,seq.name) %>%
  filter(str_detect(seq.name,"EVMZ")) %>%
  distinct() %>%
  
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EV")) %>%
      distinct(seq.name) %>%
      separate(seq.name, into = c("seq.name", "isomer"), sep = "-") 
  ) %>% 
  mutate(seq_name = str_c(seq.name, isomer, sep = "-")) %>%
  distinct(seq_name) %>%
  rename(seq.name = seq_name) %>%
  left_join(
    EV_MAB_prot_fasta %>%
      filter(str_detect(seq.name,"EV")) %>%
      select(-genome)
  ) %>%
  # rename(seq.name = seq_name) %>%
  # mutate(seq.text = str_replace_na(seq.text, "NA")) %>%
  # filter(seq.text == "NA")
  mutate(seq.name = paste0(">", seq.name)) %>%
  select(seq.name,seq.text) %>%
  pivot_longer(cols = c(seq.name, seq.text)) %>%
  select(-name) %>%
  distinct() %>%
  write.table(paste0(result_dir, "EV_specific_DNA_intergration_transposable.prot.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)
