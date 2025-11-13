# set working directory 
setwd("D:/Sadik/EV_MAB_coverage/")


# load library 
library(tidyverse)

library(karyoploteR)
library(GenomicRanges)

# Load data 

# load("bedadeti_MAB_mapping_coverage.RData")
Musa_acuminata.gene_complement.cov <-
read.delim("Musa_acuminata.gene_complement.cov.bed", header = F) %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V7) %>%
  select(chr, start, end, coverage) %>%
  mutate(coverage = coverage*100) %>%
  filter(coverage < 1)

# head(bedadeti_musa_ba_align_cov)

# format input data and generate grange 

## bedadeti

## summary of log range 

Musa_acuminata.gene_complement.cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  summarise(max = max(coverage),
            min = min(coverage),
            mean = mean(coverage))


# generate grange 

mazia_Musa_acuminata.gene_complement.cov_grange <-
  Musa_acuminata.gene_complement.cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

mazia_Musa_acuminata.gene_complement.cov_grange %>%
  head()

# remover unncessary files 
# rm(bedadeti_musa_ac_align_cov)

#  Bedadeti 

bedadeti_Musa_acuminata.gene_complement.cov <-
  read.delim("Bedadeti_Musa_acuminata.gene_complement.cov.bed", header = F) %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                coverage = V7) %>%
  select(chr, start, end, coverage) %>%
  mutate(coverage = coverage*100) %>%
  filter(coverage < 1)



# Summary of log range 
bedadeti_Musa_acuminata.gene_complement.cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  summarise(max = max(coverage),
            min = min(coverage),
            mean = mean(coverage))

# max  min     mean
# 1 11.94 0.04 4.946374

# Generate grange

bedadeti_Musa_acuminata.gene_complement.cov_grange <-
  bedadeti_Musa_acuminata.gene_complement.cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  # filter(coverage > 1) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


bedadeti_Musa_acuminata.gene_complement.cov_grange %>%
  tail



# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges

MA_genome_bed <-
  read.delim("Musa_acuminata.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name)

centro_density <-
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
  bed_to_granges ()


kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_Musa_acuminata.gene_complement.cov_grange,
       y1=mazia_Musa_acuminata.gene_complement.cov_grange$id,
       ymax=max(mazia_Musa_acuminata.gene_complement.cov_grange$id)*.2,
       col="#006400", r0=0.7, r1=.3, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_Musa_acuminata.gene_complement.cov_grange,
       y1=bedadeti_Musa_acuminata.gene_complement.cov_grange$id,
       ymax=max(bedadeti_Musa_acuminata.gene_complement.cov_grange$id*.2),
       col="#0000FF", r0=0.1, r1=.3,border=NA)


# mostdepth coverage


mazia_musa_ac_align_cov <- 
  read.delim("/Sadik/Ensete_Musa_coverage/mosdepth_coverage/mosdepth_coverage/mosdepth_coverage_mazia_musa_ac_10k_non_overlapping.bed.regions.bed.gz",
             header = F ) %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3, 
                coverage =V4) %>%
  mutate(chr = str_remove(chr,"B"))

bedadeti_musa_ac_align_cov <- 
  read.delim("/Sadik/Ensete_Musa_coverage/mosdepth_coverage/mosdepth_coverage/mosdepth_coverage_bedadeti_musa_ac_10k_non_overlapping.bed.regions.bed.gz",
             header = F ) %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3, 
                coverage =V4) %>%
  mutate(chr = str_remove(chr,"B"))

bedadeti_musa_ac_align_cov_grange <-
  bedadeti_musa_ac_align_cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

bedadeti_musa_ac_align_cov_grange %>%
  head()

rm(bedadeti_musa_ac_align_cov)


# Mazia vs Musa/Ensete_gl

# load("mazia_EG_MABS_alignment_coverage_dedup.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1)

# rm(mazia_musa_ba_align_cov,mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

# mazia_musa_ac_align_cov %>%
#   head()

# Generate grange

mazia_musa_ac_align_cov_filtered_grange <-
  mazia_musa_ac_align_cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 1) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


mazia_musa_ac_align_cov_filtered_grange %>%
  tail

mazia_musa_ac_align_cov %>% head

rm(mazia_musa_ac_align_cov)

# Generate karyotype plot

library(karyoploteR)
library(GenomicRanges)
# Genomic ranges

MA_genome_bed <-
  read.delim("Musa_acuminata.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name)

centro_density <-
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
  bed_to_granges ()


kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_musa_ac_align_cov_filtered_grange,
       y1=mazia_musa_ac_align_cov_filtered_grange$id,
       ymax=max(mazia_musa_ac_align_cov_filtered_grange$id)/2,
       col="#006400", r0=0.6, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.4,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_musa_ac_align_cov_grange,
       y1=bedadeti_musa_ac_align_cov_grange$id,
       ymax=max(bedadeti_musa_ac_align_cov_grange$id)/2,
       col="#0000FF", r0=0.41, r1=0,border=NA)



