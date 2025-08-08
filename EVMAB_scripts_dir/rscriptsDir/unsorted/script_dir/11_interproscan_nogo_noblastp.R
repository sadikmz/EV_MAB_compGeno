# read.delim(paste0(DATA_DIR,"/EV_MAB.nogo_noreads.atha_osa.diamond.blastp.e05.out.gz"))%>% 
  
  list.files(paste0(DATA_DIR,"interproscan_tsv/"), pattern = "*.tsv")

  
  EV_MAB_interprocan <- c()
  
  for(i in list.files(paste0(DATA_DIR,"interproscan_tsv/"), pattern = "*.tsv")){
    
    BASENAME <- i %>% 
      str_remove(".interprocan.01.tsv") 

    
    EV_MAB_interprocan <- 
      read.delim(paste0(DATA_DIR,"interproscan_tsv/",i), header = F) %>% 
      # mutate(genome = BASENAME) %>%
      rbind(EV_MAB_interprocan)
  }
  
  EV_MAB_interprocan <- 
  EV_MAB_interprocan %>%
    rename(seq.name = V1) %>%
    separate(seq.name, into = c("seq.name","isomer"), sep="-")

  