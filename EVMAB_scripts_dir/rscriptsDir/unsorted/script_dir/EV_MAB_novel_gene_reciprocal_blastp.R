


read.delim(paste0("Z:/sadik/evmab_dir/","EV_MAB.novel_genes.self_bp.blastp.e05.usensitive.tsv.gz"), header = F) %>%
  rename(
    qseqid = "V1",
    qlen  = "V2",
    sallseqid = "V3",
    slen = "V4",
    qstart = "V5",
    qend = "V6", 
    sstart = "V7",
    send = "V8",
    evalue = "V9",
    bitscore = "V10",
    length = "V11",
    pident = "V12",
    nident = "V13",
    mismatch = "V14",
    gapopen = "V15",
    qcovhsp = "V16",
    scovhsp = "V17",
  ) %>% 
  mutate(self_hit = str_detect(qseqid,sallseqid)) %>%
  filter(self_hit != "TRUE") %>%
  # filter(str_detect(qseqid,"EVMZ")) %>% 
  filter(pident > 80, qcovhsp > 80) %>% 
  distinct(qseqid,sallseqid) %>%
  write.table(paste0("Z:/sadik/evmab_dir/","EV_MAB.novel_genes.self_bp.blastp.e05.8080.txt"), 
              col.names =  T,
              row.names = F,
              sep = '\t',
              quote = F) 
  
  # rbind(
  #   
  #   read.delim(paste0("Z:/sadik/evmab_dir/","EV_MAB.novel_genes.self_bp.blastp.e05.usensitive.tsv.gz"), header = F) %>%
  #     rename(
  #       qseqid = "V1",
  #       qlen  = "V2",
  #       sallseqid = "V3",
  #       slen = "V4",
  #       qstart = "V5",
  #       qend = "V6", 
  #       sstart = "V7",
  #       send = "V8",
  #       evalue = "V9",
  #       bitscore = "V10",
  #       length = "V11",
  #       pident = "V12",
  #       nident = "V13",
  #       mismatch = "V14",
  #       gapopen = "V15",
  #       qcovhsp = "V16",
  #       scovhsp = "V17",
  #     ) %>% 
  #     mutate(self_hit = str_detect(qseqid,sallseqid)) %>%
  #     filter(self_hit != "TRUE") %>%
  #     filter(pident > 90, qcovhsp > 80) %>%
  #     distinct(sallseqid) %>%
  #     rename(qseqid =sallseqid )
  # ) %>%
  # distinct(qseqid) %>%
  #   filter(str_detect(qseqid,"EVMZ"))
  # 
