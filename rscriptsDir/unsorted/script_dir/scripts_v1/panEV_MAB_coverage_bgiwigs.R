# set working directory 
setwd("D:/Sadik/EV_MAB_coverage//")


# load library 
library(tidyverse)

library(karyoploteR)
library(GenomicRanges)

# Load data 

load("PanEV_MAB_mapping_coverage.RData")
# Musa_acuminata.gene_complement.cov <-
#   read.delim("Musa_acuminata.gene_complement.cov.bed", header = F) %>%
#   dplyr::rename(chr = V1,
#                 start = V2,
#                 end = V3,
#                 coverage = V7) %>%
#   select(chr, start, end, coverage) %>%
#   mutate(coverage = coverage*100) %>%
#   filter(coverage < 1)

# head(bedadeti_musa_ba_align_cov)

# format input data and generate grange 

## bedadeti

## summary of log range 

PanEV_musa_ac_align_cov %>% head()
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

panEV_Musa_acuminata.gene_complement.cov_grange <-
  PanEV_musa_ac_align_cov %>% 
  # filter(chr == "chr01") %>%
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
  # filter(coverage > 0) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()

# 
# #  for chr1
# 
# 
# panEV_Musa_acuminata.gene_complement.cov_grange <-
#   PanEV_musa_ac_align_cov %>% 
#   filter(chr == "chr01") %>%
#   # dplyr::rename(y = coverage) %>%
#   filter(str_detect(chr,"^chr")) %>%
#   # filter(coverage > 0) %>%
#   # mutate(
#   #   start = as.numeric(start),
#   #   coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
#   #                        TRUE ~ as.numeric(round(log(coverage),2)))) %>%
#   bed_to_granges ()


panEV_Musa_acuminata.gene_complement.cov_grange %>%
  head()





# Summary of log range 
PanEV_musa_ba_align_cov %>%
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

panEV_Musa_balbisiana.gene_complement.cov_grange_grange <-
  PanEV_musa_ba_align_cov %>% 
  filter(str_detect(chr,"^Bchr")) %>%
    mutate( chr = str_remove(chr,"B")) %>%
  # filter(chr == "chr01") %>%
  # filter(coverage > 1) %>%
  mutate(
    start = as.numeric(start),
    coverage = case_when(coverage <= 1 ~ round(log10(1.1),2),
                         TRUE ~ as.numeric(round(log(coverage),2)))) %>%
  bed_to_granges ()


panEV_Musa_balbisiana.gene_complement.cov_grange_grange %>%
  tail



# Centromere and karyotype range


centro_density <-
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
    # filter(chr == "chr01") %>%
    bed_to_granges ()



centro_density_MA <-
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
  bed_to_granges ()


centro_density_MB <- 
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^B")) %>%
  mutate(chr = str_replace(chr, "B","chr")) %>% 
  bed_to_granges () 




# Read bigwigs and reformat the data

bigwig_panEV_MA = "panEV_Musa_acuminata.coverage.bw"
bigwig_panEV_MB = "panEV_Musa_balbisiana.coverage.bw"


bigwig_panEV_MA_grange <-
  read_bigwig(bigwig_panEV_MA) %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(seqnames,"chr")) %>%
  as_granges()



bigwig_panEV_MB_grange <- 
  read_bigwig(bigwig_panEV_MB) %>%
  as.data.frame() %>%
  filter(str_detect(seqnames,"Bchr")) %>%
  mutate(seqnames = str_remove(seqnames,"B")) %>%
  dplyr::filter(str_detect(seqnames,"chr")) %>%
  as_granges()


#  Generate karyotype plot 


kp <- plotKaryotype(genome=MAB_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, panEV_Musa_acuminata.gene_complement.cov_grange,
       y1=panEV_Musa_acuminata.gene_complement.cov_grange$id,
       ymax=max(panEV_Musa_acuminata.gene_complement.cov_grange$id)/15,
       col="#006400", r0=0.55, r1=1, border="#006400")
       # col="#006400", r0=0.7, r1=.3, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="black", r0=0.42, r1=1,border=NA)

kpBars(kp, panEV_Musa_balbisiana.gene_complement.cov_grange_grange,
       y1=panEV_Musa_balbisiana.gene_complement.cov_grange_grange$id,
       ymax=max(panEV_Musa_balbisiana.gene_complement.cov_grange_grange$id/15),
       col="#0000FF", r0=0.1, r1=.3,border="#0000FF")




#########




bigwig_panEV_MA <-
  read_bigwig(bigwig_panEV_MA) %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(seqnames,"chr")) %>%
  mutate(#seqnames = str_remove(seqnames,"B"),
    score = log10(score),
    score = case_when(score=="-Inf" ~ 0,
                      TRUE ~ score)) %>%
  as_granges()


bigwig_panEV_MB <-
  read_bigwig(bigwig_panEV_MB) %>%
  as.data.frame() %>%
  filter(str_detect(seqnames,"Bchr")) %>%
  mutate(seqnames = str_remove(seqnames,"B"),
         score = log10(score),
         score = case_when(score=="-Inf" ~ 0,
                           TRUE ~ score)
         ) %>%
  dplyr::filter(str_detect(seqnames,"chr")) %>%
  as_granges()

MAB_genome_bed <-
  MAB_genome_bed %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(seqnames,"chr")) %>%
  as_granges()

ymax_uppper= ceiling(max(abs(range(bigwig_panEV_MA$score))))
ymin_uppper= -ymax_uppper

ymax_lower= ceiling(max(abs(range(bigwig_panEV_MB$score))))
ymin_lower= -ymax_lower


## Open device 

png("panEV_MAB_coverage.bigiwigs.v10.png",width = 650, height = 750)

kp <- plotKaryotype(genome=MAB_genome_bed,  ideogram.plotter = NULL, plot.type=2 )

# PanEV reads coverage against MA
kp <- kpBars(kp, bigwig_panEV_MA, 
             y1=bigwig_panEV_MA$score, 
             ymax = max(bigwig_panEV_MA$score)/4,
             col="green", r0=-0.1, r1=0.4,  border=NA )

kpAxis(kp, ymax = ymax_uppper, ymin = ymin_uppper, r1 = 0.8, cex = 0.8 )

# MA centromeric repeats 

kpBars(kp, centro_density_MA,
       y1=centro_density_MA$id,
       ymax=max(centro_density_MA$id)/1.5,
       col="red", r0=-0.3, r1=-0.1,border=NA)

# MB centromeric repeats
kpBars(kp, centro_density_MB,
       y1=centro_density_MB$id,
       ymax=max(centro_density_MB$id)/1.5,
       col="black", r0=-0.3, r1=-0.5,border=NA)


kp <- kpBars(kp, bigwig_panEV_MB, 
             y1=bigwig_panEV_MB$score, 
             ymax = max(bigwig_panEV_MB$score)/4,
             col="#0000FF", r0=-0.48, r1=-0.9, border=NA )
# kpAxis(kp, ymax = ymax_uppper, ymin = ymin_uppper, r0  = 0.8, cex = 1 )


dev.off()






