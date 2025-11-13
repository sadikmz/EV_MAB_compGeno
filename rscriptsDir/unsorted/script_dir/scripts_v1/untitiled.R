

load("bedadeti_EG_MABS_alignment_coverage.RData")


bedadti_musa_ac_align_cov <- mazia_musa_ac_align_cov
# bedadeti_musa_ac_align_cov <- mazia_musa_ba_align_cov
# bedadeti_musa_sc_align_cov <- mazia_musa_sc_align_cov
# bedadti_musa_ac_align_cov <- mazia_ensete_gl_align_cov

rm (mazia_musa_ac_align_cov,mazia_musa_ba_align_cov,mazia_musa_sc_align_cov, mazia_ensete_gl_align_cov)

bedadeti_musa_ac_align_cov %>%
  tail()

## format input data and generate grange 

# bedadeti 

bedadti_musa_ac_align_cov_grange <-
  bedadeti_musa_ac_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^Bchr")) %>%
  filter(coverage > 1) %>%
  mutate(coverage = round(log10(coverage),1),
         chr = str_replace(chr, "Bchr","chr")) %>%
  bed_to_granges ()

bedadti_musa_ac_align_cov_grange %>%
  head()

rm(bedadeti_musa_ac_align_cov)

# Mazia

load("mazia_EG_MABS_alignment_coverage.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1) 

rm(mazia_ensete_gl_align_cov,mazia_musa_ac_align_cov,mazia_musa_sc_align_cov)

# mazia_musa_ba_align_cov %>%
#   head()

# Generate grange 

mazia_musa_ba_align_cov_filtered_grange <-
  mazia_musa_ba_align_cov %>% 
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^Bchr")) %>%
  filter(coverage > 1) %>%
  mutate(coverage = round(log10(coverage),1),
         chr = str_replace(chr, "Bchr","chr")) %>%
  bed_to_granges ()



mazia_musa_ba_align_cov_filtered_grange %>%
  tail

rm (mazia_musa_ba_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges 

MB_genome_bed <- read.delim("Musa_balbisiana.genome.bed") %>%
  filter(str_detect(seq.name,"Bchr")) %>%
  dplyr::rename(chr = seq.name) %>%
  mutate(chr = str_replace(chr, "Bchr","chr"))


MB_genome_bed 

