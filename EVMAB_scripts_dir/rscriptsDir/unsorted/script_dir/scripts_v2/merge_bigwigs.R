
list.files(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/"),
           pattern = ".mapped.coverage.bw")

reads_coverage_bigwig = list()


# reads data
for (i in list.files(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/"),
                     pattern = ".mapped.coverage.bw")) {
  BASENAME <- basename(i) %>%
    # str_remove(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/", i)) %>% 
                 str_remove("_Musa_balbisiana.mapped.coverage.bw")
  
  read_bw <- 
    plyranges::read_bigwig(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME)
  
  
  reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}

# Unlist

reads_coverage_bigwig =
  do.call(rbind, reads_coverage_bigwig) 
# remove row.names()

row.names(reads_coverage_bigwig) <- NULL
rm(read_cov)


#qtrimmed 

reads_coverage_bigwig = list()


# reads data
for (i in list.files(paste0(DATA_DIR,"musa_balbisiana/"),
                     pattern = ".qfiltered.coverage.bw")) {
  BASENAME <- basename(i) %>%
    # str_remove(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/", i)) %>% 
    str_remove(".qfiltered.coverage.bw")
  
  read_bw <- plyranges::read_bigwig(paste0(DATA_DIR,"musa_balbisiana/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME)
  
  
  reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}

# Unlist

reads_coverage_bigwig =
  do.call(rbind, reads_coverage_bigwig) 
# remove row.names()

row.names(reads_coverage_bigwig) <- NULL
rm(read_cov)

####
reads_coverage_bigwig = c()


# reads data
for (i in list.files(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/"),
                     pattern = ".mapped.coverage.bw")) {
  BASENAME <- basename(i) %>%
    # str_remove(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/", i)) %>% 
    str_remove("_Musa_balbisiana.mapped.coverage.bw")
  
  reads_coverage_bigwig <- 
    plyranges::read_bigwig(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME) %>%
    rbind(reads_coverage_bigwig)
  
  
  # reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}

# Unlist

# reads_coverage_bigwig =
#   do.call(rbind, reads_coverage_bigwig) 
# remove row.names()

row.names(reads_coverage_bigwig) <- NULL


reads_coverage_bigwig %>%
  distinct(reads_source)
### qfiltered

qfiltered_reads_coverage_bigwig = c()


# reads data
for (i in c("Musa_balbisiana_arkiya_dawro.qfiltered.coverage.bw", "Musa_balbisiana_arkiya_wolaita.qfiltered.coverage.bw", 
            "Musa_balbisiana_astara.qfiltered.coverage.bw", "Musa_balbisiana_bedadeti.qfiltered.coverage.bw", 
            "Musa_balbisiana_buffero.qfiltered.coverage.bw", "Musa_balbisiana_derea.qfiltered.coverage.bw", 
            "Musa_balbisiana_epo.qfiltered.coverage.bw", "Musa_balbisiana_erpha.qfiltered.coverage.bw" 
)) {
  BASENAME <- basename(i) %>%
    str_remove("Musa_balbisiana_") %>%
    str_remove(".qfiltered.coverage.bw")
  
  qfiltered_reads_coverage_bigwig <- 
    plyranges::read_bigwig(paste0(DATA_DIR,"musa_balbisiana/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME) %>%
    rbind(qfiltered_reads_coverage_bigwig)
  
  
  # reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}
qfiltered_reads_coverage_bigwig %>%
  distinct(reads_source)

qfiltered_reads_coverage_bigwig_v1 <-c ()

for (i in c( 
  "Musa_balbisiana_nechuwe.qfiltered.coverage.bw", "Musa_balbisiana_nobo.qfiltered.coverage.bw", 
  "Musa_balbisiana_onjamo.qfiltered.coverage.bw", "Musa_balbisiana_siyuti.qfiltered.coverage.bw", 
  "Musa_balbisiana_yako.qfiltered.coverage.bw", "Musa_balbisiana_yanbule.qfiltered.coverage.bw"
)) {
  BASENAME <- basename(i) %>%
    str_remove("Musa_balbisiana_") %>%
    str_remove(".qfiltered.coverage.bw")
  
  qfiltered_reads_coverage_bigwig_v1 <- 
    plyranges::read_bigwig(paste0(DATA_DIR,"musa_balbisiana/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME) %>%
    rbind(qfiltered_reads_coverage_bigwig_v1)
  
  
  # reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}

qfiltered_reads_coverage_bigwig_v2 <- c()


for (i in c(            "Musa_balbisiana_erpha13.qfiltered.coverage.bw", "Musa_balbisiana_erpha20.qfiltered.coverage.bw", 
                        "Musa_balbisiana_lochinge_wolaita.qfiltered.coverage.bw", "Musa_balbisiana_mazia_dawro.qfiltered.coverage.bw",
            "Musa_balbisiana_mazia_wolaita.qfiltered.coverage.bw", "Musa_balbisiana_maziastlfr.qfiltered.coverage.bw")) {
  BASENAME <- basename(i) %>%
    str_remove("Musa_balbisiana_") %>%
    str_remove(".qfiltered.coverage.bw")
  
  qfiltered_reads_coverage_bigwig_v2 <- 
    plyranges::read_bigwig(paste0(DATA_DIR,"musa_balbisiana/", i)) %>%
    as.data.frame() %>%
    mutate(reads_source = BASENAME) %>%
    rbind(qfiltered_reads_coverage_bigwig_v2)
  
  
  # reads_coverage_bigwig[[i]] <- read_bw # add it to your list
  
}

# Unlist

# reads_coverage_bigwig =
#   do.call(rbind, reads_coverage_bigwig) 
# remove row.names()

row.names(qfiltered_reads_coverage_bigwig_v1) <- NULL

reads_coverage_bigwig %>%
  distinct(reads_source) %>%
  arrange(reads_source)

# qfiltered_reads_coverage_bigwig <- 
# qfiltered_reads_coverage_bigwig %>%
#   filter(reads_source!="epo",
#          reads_source!="erpha") %>%
#   distinct(reads_source) %>%
#   arrange(reads_source)

# reads_coverage_bigwig <- 
# reads_coverage_bigwig %>%
qfiltered_reads_coverage_bigwig_v3 <-
  bind_rows(
    qfiltered_reads_coverage_bigwig_v1 %>%
        filter(reads_source=="yako"|
               reads_source=="yanbule"),

    qfiltered_reads_coverage_bigwig_v2 %>%
      filter(reads_source=="maziastlfr")) %>%
  distinct() 

rm(qfiltered_reads_coverage_bigwig_v1,
   qfiltered_reads_coverage_bigwig_v2)

reads_coverage_bigwig_v1 <-
reads_coverage_bigwig %>%
bind_rows(qfiltered_reads_coverage_bigwig_v3) #%>%
  # distinct()

rm(qfiltered_reads_coverage_bigwig_v3)
rm(reads_coverage_bigwig)

reads_coverage_bigwig <- 
reads_coverage_bigwig_v1 %>%
  distinct()


rm(reads_coverage_bigwig_v1)


reads_coverage_bigwig %>%
  distinct(reads_source)

reads_coverage_bigwig  %>%
  filter(reads_source!="erpha20",
           reads_source !="erpha13") %>%
  select(-reads_source) %>%
  group_by(seqnames,start,end,strand) %>%
  summarise(width = round(mean(width),1),
            score = round(mean(score),1)) %>%
  select(seqnames, start, end, width, strand, score) %>%
  head()

rm(qfiltered_reads_coverage_bigwig)

reads_coverage_bigwig %>%
  head()
write_bigwig(x, file)


# 
# write_bigwig(reads_coverage_bigwig %>%
#                dplyr::select(-reads_source),
#              paste0(DATA_DIR,"panEV_MB_coverage_bigwig.bw"))


## read bedgraph and output average values

combined.bdgraph <- 
read.delim(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/domesticated/combined.bdgraph"), header = F)

combined.bdgraph %>% 
  group_by(V1,V2,V4) %>%
  summarise(V4 = mean(V4)) %>%
  write.table(paste0(DATA_DIR,"panEV_MAB_updated/panEV_MB/domesticated/combined.mean.bedgraph"),
              col.names = F,
              row.names = F,
              quote = F,
              sep= "\t"
              )


read.delim("combined.bdgraph", header = F) %>%
  group_by(V1,V2,V4) %>%
  summarise(V4 = mean(V4)) %>%
  write.table("combined.mean.bedgraph",
              col.names = F,
              row.names = F,
              quote = F,
              sep= "\t"
  )
  # head()






