setwd("D:/Sadik/Ensete_Musa_coverage/")

# Bedadeti vs Musa/Ensete_gl
library(tidyverse)
# load("bedadeti_EG_MABS_alignment_coverage_dedup.RData")

load("EV_EG_MABS_mapping_coverage_R_objects/mazia_MB_alignment_coverage_dedup.RData")
load("EV_EG_MABS_mapping_coverage_R_objects/bedadeti_MB_alignment_coverage_dedup.RData")


# bedadeti_musa_ac_align_cov <- mazia_musa_ac_align_cov
# bedadeti_musa_ba_align_cov <- mazia_musa_ba_align_cov
# bedadeti_musa_sc_align_cov <- mazia_musa_sc_align_cov
# bedadeti_musa_ac_align_cov <- mazia_ensete_gl_align_cov

# rm (mazia_musa_ba_align_cov,mazia_musa_ac_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

bedadeti_musa_ba_align_cov %>%
  tail()


## format input data and generate grange 

# bedadeti 

bedadeti_musa_ba_align_cov_grange <-
  bedadeti_musa_ba_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^Bchr")) %>%
  filter(coverage > 0)  %>%
  mutate(
    chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

bedadeti_musa_ba_align_cov_grange %>%
  head()

head(bedadeti_musa_ba_align_cov)

rm(bedadeti_musa_ba_align_cov)


# Mazia 

# load("mazia_EG_MABS_alignment_coverage_dedup.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1) 

# rm(mazia_musa_ac_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

mazia_musa_ba_align_cov %>%
  head()

# Generate grange 

mazia_musa_ba_align_cov_filtered_grange <-
  mazia_musa_ba_align_cov %>% 
    # head() %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^Bchr")) %>%
  filter(coverage > 0)  %>%
  mutate(
    chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


mazia_musa_ba_align_cov_filtered_grange %>%
  tail

mazia_musa_ba_align_cov %>% head

rm(mazia_musa_ba_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges 

MB_genome_bed <- 
  read.delim("Musa_balbisiana.genome.bed") %>%
  filter(str_detect(seq.name,"Bchr")) %>%
  dplyr::rename(chr = seq.name) %>%
  mutate(chr = str_replace(chr, "Bchr","chr") )


centro_density <- 
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^B")) %>%
  mutate(chr = str_replace(chr, "B","chr")) %>% 
  bed_to_granges () 


kp <- plotKaryotype(genome=MB_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_musa_ba_align_cov_filtered_grange,
       y1=mazia_musa_ba_align_cov_filtered_grange$id,
       ymax=max(mazia_musa_ba_align_cov_filtered_grange$id)/2,
       col="#006400", r0=0.55, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_musa_ba_align_cov_grange,
       y1=bedadeti_musa_ba_align_cov_grange$id,
       ymax=max(bedadeti_musa_ba_align_cov_grange$id)/3.2,
       col="#0000FF", r0=0.35, r1=0,border=NA)

# using mosdepth 

setwd("D:/Sadik/Ensete_Musa_coverage/")

# Bedadeti vs Musa/Ensete_gl
library(tidyverse)
# load("bedadeti_EG_MABS_alignment_coverage_dedup.RData")

mazia_musa_ba_align_cov <- 
  read.delim("/Sadik/Ensete_Musa_coverage/mosdepth_coverage/mosdepth_coverage/mosdepth_coverage_mazia_musa_ba_10k_non_overlapping.bed.regions.bed.gz",
           header = F ) %>%
  dplyr::rename(chr = V1,
         start = V2,
         end = V3, 
         coverage =V4) %>%
  mutate(chr = str_remove(chr,"B"))

bedadeti_musa_ba_align_cov <- 
  read.delim("/Sadik/Ensete_Musa_coverage/mosdepth_coverage/mosdepth_coverage/mosdepth_coverage_bedadeti_musa_ba_10k_non_overlapping.bed.regions.bed.gz",
             header = F ) %>%
  dplyr::rename(chr = V1,
         start = V2,
         end = V3, 
         coverage =V4) %>%
  mutate(chr = str_remove(chr,"B"))


bedadeti_musa_ba_align_cov %>%
  tail()


## format input data and generate grange 

# bedadeti 

bedadeti_musa_ba_align_cov_grange <-
  bedadeti_musa_ba_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0)  %>%
  mutate(
    chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

bedadeti_musa_ba_align_cov_grange %>%
  head()

rm(bedadeti_musa_ba_align_cov)


# Mazia 

# load("mazia_EG_MABS_alignment_coverage_dedup.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1) 

# rm(mazia_musa_ac_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

mazia_musa_ba_align_cov %>%
  head()

# Generate grange 

mazia_musa_ba_align_cov_filtered_grange <-
  mazia_musa_ba_align_cov %>% 
  # head() %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0)  %>%
  mutate(
    chr = str_replace(chr, "Bchr", "chr"),
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


mazia_musa_ba_align_cov_filtered_grange %>%
  tail

mazia_musa_ba_align_cov %>% head

rm(mazia_musa_ba_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges 

MB_genome_bed <- 
  read.delim("Musa_balbisiana.genome.bed") %>%
  filter(str_detect(seq.name,"Bchr")) %>%
  dplyr::rename(chr = seq.name) %>%
  mutate(chr = str_replace(chr, "Bchr","chr") )


centro_density <- 
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^B")) %>%
  mutate(chr = str_replace(chr, "B","chr")) %>% 
  bed_to_granges () 


kp <- plotKaryotype(genome=MB_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_musa_ba_align_cov_filtered_grange,
       y1=mazia_musa_ba_align_cov_filtered_grange$id,
       ymax=max(mazia_musa_ba_align_cov_filtered_grange$id)/1.75,
       col="#006400", r0=0.6, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.4,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_musa_ba_align_cov_grange,
       y1=bedadeti_musa_ba_align_cov_grange$id,
       ymax=max(bedadeti_musa_ba_align_cov_grange$id)/1.75,
       col="#0000FF", r0=0.41, r1=0,border=NA)


save.image(file = "EV_MAB_coverage.RData")
load("EV_MAB_coverage.RData")


