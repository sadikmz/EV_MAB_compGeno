


bigwig_panEV_MA = "bigwigs/panEV_Musa_acuminata.coverage.bw"
bigwig_panEV_MB = "panEV_Musa_balbisiana.coverage.bw"


bigwig_panEV_MA <-
  read_bigwig(bigwig_panEV_MA) %>%
  as.data.frame() #%>%
  # dplyr::filter(str_detect(seqnames,"chr")) %>%
  # as_granges()


bigwig_panEV_MA %>%
  filter(seqnames == "chr04") %>%
  # head()
  filter(start > 4457583 ) %>%
  filter(end < 5950000) %>%
  tail()



# 
# bigwig_panEV_MB_grange <- 
#   read_bigwig(bigwig_panEV_MB) %>%
#   as.data.frame() %>%
#   filter(str_detect(seqnames,"Bchr")) %>%
#   mutate(seqnames = str_remove(seqnames,"B")) %>%
#   dplyr::filter(str_detect(seqnames,"chr")) %>%
#   as_granges()



panEV_MA_genomecov <-
  read.delim("bigwigs/panEV_Musa_acuminata.genomecov.bed", header = F)

panEV_MA_genomecov %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >0 ) %>%
  summarise(tot_length = sum(length))


panEV_MA_genomecov_gene <-
  read.delim("bigwigs/panEV_Musa_acuminata.genomecov.gene.bed", header = F)

panEV_MA_genomecov_gene %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >0 ) %>%
  summarise(tot_length = sum(length))




panEV_MB_genomecov <-
  read.delim("bigwigs/panEV_Musa_balbisiana.genomecov.v1.bed", header = F)

panEV_MB_genomecov %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >0 ) %>%
  summarise(tot_length = sum(length))



panEV_MB_genomecov_gene <-
  read.delim("bigwigs/panEV_Musa_balbisiana.genomecov.gene.v1.bed", header = F)

panEV_MB_genomecov_gene %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  # filter(coverage >0 ) %>%
  summarise(tot_length = sum(length))


panAA_mazia_genomecov <-
  read.delim("bigwigs/panAA_EV_mazia.allMapped.coverage.bed", header = F)

panAA_mazia_genomecov %>% 
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >=1 ) %>%
  summarise(tot_length = sum(length))



panAA_bedadeti_genomecov <-
  read.delim("bigwigs/panAA_EV_bedadeti.allMapped.coverage.bed", header = F)

panAA_bedadeti_genomecov %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >0 ) %>%
  summarise(tot_length = sum(length))




panBB_mazia_genomecov <-
  read.delim("bigwigs/panBB_EV_mazia.allMapped.coverage.bed", header = F)

panBB_mazia_genomecov %>% 
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >=1 ) %>%
  summarise(tot_length = sum(length))



panBB_bedadeti_genomecov <-
  read.delim("bigwigs/panBB_EV_bedadeti.allMapped.coverage.bed", header = F)

panBB_bedadeti_genomecov %>% 
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V4) %>%
  select(chr, start, end, coverage) %>%
  mutate(length = end - start ) %>%
  filter(coverage >=1 ) %>%
  summarise(tot_length = sum(length))




