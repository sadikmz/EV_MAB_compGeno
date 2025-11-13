
# Bedadeti vs Musa/Ensete_gl
library(tidyverse)
# load("bedadeti_EG_MABS_alignment_coverage_dedup.RData")
setwd("D:/Sadik/Ensete_Musa_coverage/")

load("EV_EG_MABS_mapping_coverage_R_objects/mazia_EG_alignment_coverage_dedup.RData")
load("EV_EG_MABS_mapping_coverage_R_objects/bedadeti_EG_alignment_coverage_dedup.RData")



# bedadti_musa_ac_align_cov <- mazia_musa_ac_align_cov
# bedadeti_musa_ba_align_cov <- mazia_musa_ba_align_cov
# bedadeti_musa_sc_align_cov <- mazia_musa_sc_align_cov
# bedadeti_ensete_gl_align_cov <- mazia_ensete_gl_align_cov


# rm (mazia_musa_ac_align_cov,mazia_musa_sc_align_cov, mazia_musa_ba_align_cov, mazia_ensete_gl_align_cov,bedadti_musa_ac_align_cov)

# bedadeti_ensete_gl_align_cov %>%
#   tail()

## format input data and generate grange 

# bedadeti 

bedadeti_ensete_gl_align_cov_grange <-
  bedadeti_ensete_gl_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

bedadeti_ensete_gl_align_cov_grange %>%
  head()




rm(bedadeti_ensete_gl_align_cov)


# read repeat file 

# read.delim("C:/Users/Sadik.DESKTOP-UQ756TB/Downloads/EG_TEs.bed") %>%
#   head()



# Mazia vs Musa/Ensete_gl

# load("mazia_EG_MABS_alignment_coverage_dedup.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1) 

# rm(mazia_musa_ac_align_cov,mazia_musa_ba_align_cov,mazia_musa_sc_align_cov)

# mazia_musa_ac_align_cov %>%
#   head()

# mazia_ensete_gl_align_cov %>%
#   head()
# Generate grange 

mazia_ensete_gl_align_covd_grange <-
  mazia_ensete_gl_align_cov %>% 
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  # filter(coverage > 1) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


mazia_ensete_gl_align_covd_grange %>%
  tail

mazia_ensete_gl_align_cov %>% head

rm(mazia_ensete_gl_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges 

EG_genome_bed <- 
  read.delim("Ensete_glaucum.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name)

TE_bed <-
  
  read.delim("ensete_gl.TE.bed", header= F, sep = '\t') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3
                #density = V4
                ) %>% 
  # dplyr::filter(str_detect(chr,"^A")) %>%
  # mutate(chr = str_replace(chr, "A","chr")) %>% 
  bed_to_granges () 


kp <- plotKaryotype(genome=EG_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_ensete_gl_align_covd_grange,
       y1=mazia_ensete_gl_align_covd_grange$id,
       ymax=max(mazia_ensete_gl_align_covd_grange$id)/2,
       col="#006400", r0=0.55, r1=1, border=NA)
# 
# kpBars(kp, centro_density,
#        y1=centro_density$id,
#        ymax=max(centro_density$id)/0.2,
#        col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_ensete_gl_align_cov_grange,
       y1=bedadeti_ensete_gl_align_cov_grange$id,
       ymax=max(bedadeti_ensete_gl_align_cov_grange$id)/3.2,
       col="#0000FF", r0=0.35, r1=0,border=NA)


