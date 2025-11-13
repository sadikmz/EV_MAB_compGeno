library("dplyr")
bedadeti_musa_ba_align_cov %>%
  mutate(len = end - start + 1) %>%
  filter( coverage < 1) %>%
  select(-len) %>%
  write.table("bedadeti_musa_ba_align_cov_0.bed",
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t')
  
 
  bedadeti_musa_ac_align_cov %>%
    mutate(len = end - start + 1) %>%
    filter( coverage < 1) %>%
    select(-len) %>%
  write.table("bedadeti_musa_ac_align_cov_0.bed",
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t')
  
  
  # head()



mazia_musa_ba_align_cov %>%
  mutate(len = end - start + 1) %>%
  filter( coverage < 1) %>%
  select(-len) %>%
  write.table("mazia_musa_ba_align_cov_0.bed",
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t')


mazia_musa_ac_align_cov %>%
  mutate(len = end - start + 1) %>%
  filter( coverage < 1) %>%
  select(-len) %>%
  write.table("mazia_musa_ac_align_cov_0.bed",
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t')
# head()


data_path ="D:/Sadik/EV_MAB_coverage/coverage_noreads//"


read.delim(paste0(data_path,"bedadeti_Musa_balbisiana.noreads.coverage.bed"), header=F) %>%
  filter(V3==0) %>%
  head
