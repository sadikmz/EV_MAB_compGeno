# load library 
library(plyranges)
library(tidyverse)
library(karyoploteR)
library(regioneR)
library(zoo)
library(karyoploteR)
library(GenomicRanges)
library(plyranges)

setwd("D:/Sadik/EV_MAB_coverage//")

set.seed(1234)

#Parameters
data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")



#

########


# bigwig_BD_MA.file <- "bedadeti_Musa_acuminata.mapped.sorted.makrdup.bw"
# 
# bigwig_BD_MA.bw <- read_bigwig(bigwig_BD_MA.file)
# 
# bigwig_MZ_MA.file <- "mazia_Musa_acuminata.mapped.sorted.makrdup.bw"
# 
# bigwig_MZ_MA.bw <- read_bigwig(bigwig_MZ_MA.file)
# 
# 
# 
# bigwig_panEV_MA_1k = "panEV_Musa_acuminata.coverage.bw"
# bigwig_panEV_MB_1k = "panEV_Musa_balbisiana.coverage.bw"
# 
# 
# bigwig_panEV_MB_1k.v1 <-
# read_bigwig(bigwig_panEV_MB_1k) %>%
#   as.data.frame() %>%
#   filter(str_detect(seqnames,"Bchr")) %>%
#   mutate(seqnames = str_remove(seqnames,"B")) %>%
#   as_granges() #%>%
#   
# export(bigwig_panEV_MB_1k.v1,con = "panEV_Musa_balbisiana.coverage.v1.bw", format = "bigWig")
# 
# export.bw("panEV_Musa_balbisiana.coverage.v1.bw")
#   mutate()
# 
# panEV_Musa_balbisiana.coverage.bw


# chromosome labesl 

MA_genome_bed <-
  read.delim("Musa_acuminata.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name) 
# filter(chr == "chr01")


Mb_genome_bed <-
  read.delim("Musa_balbisiana.genome.bed") %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(chr = seq.name) 
# filter(chr == "chr01")


#  Chromosome length adjusted to both MA and MB
MAB_genome_bed <-
  MA_genome_bed %>%
  left_join(Mb_genome_bed %>%
              mutate(chr= str_remove (chr,"B")) %>%
              dplyr::rename(end_MB = "end")) %>%
  mutate(end_MAB = case_when(end_MB > end ~ end_MB,
                             TRUE ~ end)) %>%
  select(-end,-end_MB) %>%
  dplyr::rename(end = end_MAB) %>%
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


  
  
  bigwig_panEV_MA = "panEV_Musa_acuminata.coverage.bw"
  bigwig_panEV_MB = "panEV_Musa_balbisiana.coverage.bw"
  

  # 
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
  
    
    # 
    
    kp <- plotKaryotype(genome=MAB_genome_bed,  ideogram.plotter = NULL, plot.type=2 )
    kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA, data.panel = 1)
    kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA, data.panel = 1)
    
    
    bigwig_panEV_MA_chr1 <-
    read_bigwig(bigwig_panEV_MA) %>%
      as.data.frame() %>%
      dplyr::filter(str_detect(seqnames,"chr01")) %>%
      mutate(#seqnames = str_remove(seqnames,"B"),
             score = log10(score),
             score = case_when(score=="-Inf" ~ 0,
                               TRUE ~ score)) %>%
      as_granges()
    
    
    bigwig_panEV_MB_chr1 <-
      read_bigwig(bigwig_panEV_MB) %>%
      as.data.frame() %>%
      filter(str_detect(seqnames,"Bchr")) %>%
      mutate(seqnames = str_remove(seqnames,"B"),
             score = log10(score),
             score = case_when(score=="-Inf" ~ 0,
              TRUE ~ score)) %>%
      dplyr::filter(str_detect(seqnames,"chr01")) %>%
      as_granges()
    
    MAB_genome_bed_chr1 <-
    MAB_genome_bed %>%
      as.data.frame() %>%
      dplyr::filter(str_detect(seqnames,"chr01")) %>%
      as_granges()
    
    ymax_uppper= ceiling(max(abs(range(bigwig_panEV_MA_chr1$score))))
    ymin_uppper= -ymax_uppper
    
    ymax_lower= ceiling(max(abs(range(bigwig_panEV_MB_chr1$score))))
    ymin_lower= -ymax_lower
    
    
    kp <- plotKaryotype(genome=MAB_genome_bed_chr1,  ideogram.plotter = NULL, plot.type=2 )
    
    # PanEV reads coverage against MA
    kp <- kpBars(kp, bigwig_panEV_MA_chr1, 
                 y1=bigwig_panEV_MA_chr1$score, 
                 ymax = max(bigwig_panEV_MA_chr1$score),
                 col="green", r0=-0.2, r1=0.3,  border=NA )
    
    kpAxis(kp, ymax = ymax_uppper, ymin = ymin_uppper, r1 = -0.2 )
    
    # MA centromeric repeats 
    
    kpBars(kp, centro_density_MA,
           y1=centro_density_MA$id,
           ymax=max(centro_density_MA$id),
           col="red", r0=-0.3, r1=-0.2,border=NA)
    
    # MB centromeric repeats
    kpBars(kp, centro_density_MB,
           y1=centro_density_MB$id,
           ymax=max(centro_density_MB$id),
           col="black", r0=-0.3, r1=-0.4,border=NA)
    
    
    kp <- kpBars(kp, bigwig_panEV_MB_chr1, 
                 y1=bigwig_panEV_MB_chr1$score, 
                 ymax = max(bigwig_panEV_MB_chr1$score),
                 col="#0000FF", r0=-0.38, r1=-0.8, border=NA )
    
    kpAxis(kp, ymax = ymax_lower, ymin = ymin_lower,r1 = -0.4 )
    

    
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
  
    kp <- plotKaryotype(genome=MAB_genome_bed_chr1,  ideogram.plotter = NULL, plot.type=2 )
    kp <- kpPlotCoverage(kp, data=bigwig_panEV_MA_chr1v, data.panel = 1)
    kp <- kpPlotCoverage(kp, data=bigwig_panEV_MA_chr1v, )
    kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA, data.panel = 1)
    kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA, data.panel = 1)
    kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA, data.panel = 1)
  
  
  
  
  
  
  bigwig_panEV_MA.v1 <-
    read_bigwig(bigwig_panEV_MA) %>%
    as.data.frame() %>%
    filter(str_detect(seqnames,"chr")) %>%
    # mutate(score = log10(score),
    #        score = case_when(score=="-Inf" ~ 0,
                             # TRUE ~ score)) %>%
    as_granges() #%>%  
  
  bigwig_panEV_MB.v1 <-
    read_bigwig(bigwig_panEV_MB) %>%
    as.data.frame() %>%
    filter(str_detect(seqnames,"Bchr")) %>%
    mutate(seqnames = str_remove(seqnames,"B")) %>%
    # mutate(score = log10(score),
    #        score = case_when(score=="-Inf" ~ 0,
    #                          TRUE ~ score)) %>%
    as_granges() #%>%
  
  
  ## Open device 
  
  png("panEV_MAB.png",width = 1500, height = 700)
  
  kp <- plotKaryotype(genome=MAB_genome_bed,  ideogram.plotter = NULL, plot.type=2 )
  
  
  kp <- kpBars(kp, bigwig_panEV_MA_1k.cv1, 
               y1=bigwig_panEV_MA_1k.cv1$score, 
               ymax = max(bigwig_panEV_MA_1k.cv1$score)/0.4,
               border="green", lwd=4, data.panel = 1, r=0.5, r1=0)
  kpBars(kp, centro_density_MA,
         y1=centro_density_MA$id,
         ymax=max(centro_density_MA$id),
         col="red", r0=-0.2, r1=0,border=NA)
  kpBars(kp, centro_density_MB,
         y1=centro_density_MB$id,
         ymax=max(centro_density_MB$id),
         col="black", r0=-0.2, r1=-0.45,border=NA)
  #
  kp <- kpBars(kp, bigwig_panEV_MB_1k.v1, 
                     y1=bigwig_panEV_MA_1k.cv1$score, 
                     ymax = max(bigwig_panEV_MA_1k.cv1$score),
                     border="blue", lwd=4, data.panel = 2, r=0, r1=3)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=kp$plot$ymin, ymax=kp$plot$ymax)
  kpAxis(kp, data.panel = 1)
  kpAxis(kp, data.panel = 2)
  kpAddBaseNumbers(kp, units = "Mb")
  
  
  ## Close teh device
  
  dev.off()
