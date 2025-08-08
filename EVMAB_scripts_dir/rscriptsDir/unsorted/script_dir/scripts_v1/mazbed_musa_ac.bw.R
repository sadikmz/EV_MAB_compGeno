# set working directory 
setwd("D:/Sadik/EV_MAB_coverage/")


# load library 
library(tidyverse)

library(karyoploteR)
library(GenomicRanges)


# test 


# if (.Platform$OS.type != "windows") {
  bigwig.file <- system.file("extdata", "BRCA.genes.hg19.bw", package = "karyoploteR")
  brca.genes.file <- system.file("extdata", "BRCA.genes.hg19.txt", package = "karyoploteR")
  brca.genes <- toGRanges(brca.genes.file)
  seqlevelsStyle(brca.genes) <- "UCSC"
  
  kp <- plotKaryotype(zoom = brca.genes[1])
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0, r1=0.2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.25, r1=0.45, border="red", lwd=2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.5, r1=0.7, ymin=0, ymax=1000, border="gold", col=NA)
  kpAxis(kp, r0=0.5, r1=0.7, ymin=0, ymax=1000)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.75, r1=0.95, ymin=0, ymax="visible.region", border="orchid", col=NA)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
# }

  
library(plyranges)
  


  bigwig_BD_MA.file <- "bedadeti_Musa_acuminata.mapped.sorted.makrdup.bw"
  
  bigwig_BD_MA.bw <- read_bigwig(bigwig_BD_MA.file)
  
  bigwig_MZ_MA.file <- "mazia_Musa_acuminata.mapped.sorted.makrdup.bw"
  
  bigwig_MZ_MA.bw <- read_bigwig(bigwig_MZ_MA.file)
  
  
  brca.genes.file <- system.file("extdata", "BRCA.genes.hg19.txt", package = "karyoploteR")
  brca.genes <- toGRanges(brca.genes.file)
  seqlevelsStyle(brca.genes) <- "UCSC"
  
  kp <- plotKaryotype(zoom = brca.genes[1])
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0, r1=0.2, )
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.25, r1=0.45, border="red", lwd=2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.5, r1=0.7, ymin=0, ymax=1000, border="gold", col=NA)
  kpAxis(kp, r0=0.5, r1=0.7, ymin=0, ymax=1000)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.75, r1=0.95, ymin=0, ymax="visible.region", border="orchid", col=NA)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
  


  
  
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
  
  
  kp <- kpPlotBigWig(kp, data=bigwig_MZ_MA.file, lwd=2, border = "#006400", data.panel = ``)
  kpBars(kp, centro_density,
         y1=centro_density$id,
         ymax=max(centro_density$id)/0.2,
         col="red", r0=0.42, r1=1,border=NA)
  kp <- kpPlotBigWig(kp, data=bigwig_BD_MA.file, r0=0.25, r1=2.5, border="#0000FF", lwd=2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.25, r1=0.45, border="red", lwd=2)
  
  # kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.5, r1=0.7, ymin=0, ymax=1000, border="gold", col=NA)
  kpAxis(kp, r0=0.5, r1=0.7, ymin=0, ymax=1000)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.75, r1=0.95, ymin=0, ymax="visible.region", border="orchid", col=NA)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
  
  
  kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
  
  
  
  # kpAddCytobandsAsLine(kp)
  kpBars(kp, bigwig_MZ_MA.bw,
         y1=bigwig_MZ_MA.bw$score,
         ymax=max(bigwig_MZ_MA.bw$score)*2,
         col="#006400", r0=0.55, r1=1, border=NA)
  
  kpBars(kp, centro_density,
         y1=centro_density$id,
         ymax=max(centro_density$id)/0.2,
         col="red", r0=0.42, r1=1,border=NA)
  
  kpBars(kp, bedadeti_musa_ac_align_cov_grange,
         y1=bedadeti_musa_ac_align_cov_grange$id,
         ymax=max(bedadeti_musa_ac_align_cov_grange$id)/3.2,
         col="#0000FF", r0=0.35, r1=0,border=NA)
  
  
  kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
  
  
  kpPlotDensity(kp, bigwig_MZ_MA.bw,
                window.size = 100,
                data.panel = "ideogram",
                col="#006400", r0=0.5, r1=1, border=NA)
  
  kpBars(kp, centro_density,
         y1=centro_density$id,
         ymax=max(centro_density$id)/0.2,
         col="red", r0=0.42, r1=1,border=NA)
  
  kpBars(kp, bedadeti_musa_ac_align_cov_grange,
         y1=bedadeti_musa_ac_align_cov_grange$id,
         ymax=max(bedadeti_musa_ac_align_cov_grange$id)/3.2,
         col="#0000FF", r0=0.35, r1=0,border=NA)
  
  
  
  
  
  # 
  
  
  
  

# read bigwig coverage track

read.

# load("bedadeti_MAB_mapping_coverage.RData")


head(bedadeti_musa_ac_align_cov)
# head(bedadeti_musa_ba_align_cov)

# format input data and generate grange 

## bedadeti

## summary of log range 

bedadeti_musa_ac_align_cov %>%
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

# remover unncessary files 
rm(bedadeti_musa_ac_align_cov)


# Mazia_MAB

# Import data 
load("mazia_MAB_mapping_coverage.RData")
# load("mazia_MAB_mapping_coverage.RData")
head(mazia_musa_ac_align_cov)
head(mazia_musa_ba_align_cov)
# rm(mazia_musa_ba_align_cov)
# rm(mazia_musa_ac_align_cov)


# Summary of log range 
mazia_musa_ac_align_cov %>%
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

mazia_musa_ac_align_cov_filtered_grange <-
  mazia_musa_ac_align_cov %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  # filter(coverage > 1) %>%
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
       col="#006400", r0=0.55, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadeti_musa_ac_align_cov_grange,
       y1=bedadeti_musa_ac_align_cov_grange$id,
       ymax=max(bedadeti_musa_ac_align_cov_grange$id)/3.2,
       col="#0000FF", r0=0.35, r1=0,border=NA)


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



