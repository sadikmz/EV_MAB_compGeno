list.files(path = paste0(DATA_DIR, "../genome_dir/"), pattern = "gff")
ma_gff <- paste0(DATA_DIR, "../genome_dir/ma.gff")

ma_gff3 <- rtracklayer::readGFF(ma_gff)
ma_gff3 %>%
  str()
# rtracklayer::export(object = "test", format = "gff3", index = TRUE)

ma_gff3 <- plyranges::read_gff3(ma_gff)

test <- ma_gff3 %>% 
  as.data.frame() %>%
  filter(type == "gene" |
           type == "mRNA") %>%
  head(n=20)


test <- ma_gff3 %>% 
  as.data.frame() %>%
  filter(type == "gene" |
           type == "mRNA") %>%
  head(n=20)

ma_gff3 %>%
  as.data.frame() %>%
  filter(str_detect(type,"mRNA")) %>%
  distinct(Note) %>%
  head()

file_gff3 <- ma_gff3 %>%
  as.data.frame()

for (i in seq_along(file_gff3$seqnames)[1:(length(seq_along(file_gff3$seqnames)) -1)]){
  
  print(i)
  
  if (file_gff3[1,7] == "gene") {
    
    # file_gff3[1,13] == file_gff3[i+1]
    
    value <- unlist(file_gff3[i+1,13])
    
    if (length(value) != 0){
      file_gff3[i,13] = unlist(file_gff3[i+1,13])
      
    }
    
    
  } else {
    print(i)
  }
}


file_gff3 %>%
  # head()
  # filter(type == "gene" |
  #          type == "mRNA") %>%
  plyranges::as_granges() %>%
  plyranges::write_gff3(paste0(DATA_DIR, "../genome_dir/ma.02.gff"))
  
  plyranges::write_gff3(paste0(DATA_DIR, "../genome_dir/ma.01.gff"), index = TRUE)

list.files(path = paste0(DATA_DIR, "../genome_dir/"),pattern = ".gff")
