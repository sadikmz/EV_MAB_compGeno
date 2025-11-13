library(tidyverse)
library(ggtext)

setwd("C:/Users/Sadik/")

load("bedadeti_EG_MABS_alignment_coverage.RData")


# bedadti_musa_ac_align_cov <- mazia_musa_ac_align_cov 
bedadeti_musa_ba_align_cov <- mazia_musa_ba_align_cov
bedadeti_musa_sc_align_cov <- mazia_musa_sc_align_cov
# bedadti_musa_ac_align_cov <- mazia_ensete_gl_align_cov

rm (mazia_musa_ac_align_cov,mazia_musa_sc_align_cov, mazia_musa_ba_align_cov, mazia_ensete_gl_align_cov)

bedadti_musa_ba_align_cov %>%
  tail()

## format input data and generate grange 

# bedadeti 

bedadti_musa_ac_align_cov_grange <-
  bedadti_musa_ac_align_cov %>% 
  # dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^chr")) %>%
  filter(coverage > 1,
         chr == "chr01" | chr=="chr02" | chr =="chr03") %>%
  mutate(coverage = round(log10(coverage),1),
         chr = str_remove(chr,"0"),
         chr = str_replace(chr, "chr","Chr")
         
         # coverage = str_replace(coverage, "-Inf","0"),
         # coverage = as.numeric(coverage)
         ) %>%
  bed_to_granges ()

bedadti_musa_ac_align_cov_grange %>%
  head()

rm(bedadti_musa_ac_align_cov)

load("mazia_EG_MABS_alignment_coverage.RData")

# mazia_musa_ac_align_cov %>%
#   filter(coverage > 1) 

# rm(mazia_ensete_gl_align_cov,mazia_musa_ba_align_cov,mazia_musa_sc_align_cov)

# mazia_ensete_gl_align_cov %>%
#   head()
# maiza Generate grange 


  mazia_musa_ac_align_cov_filtered_grange_chr1_3 <-
  mazia_musa_ac_align_cov %>% 
  # dplyr::rename(y = coverage) %>%
  filter(str_detect(chr,"^chr")) %>%
    filter(coverage > 1,
           chr == "chr01" | chr=="chr02" | chr =="chr03") %>%  
    mutate(
    coverage = round(log10(coverage),1),
         chr = str_remove(chr,"0"),
         chr = str_replace(chr, "chr","Chr"),
         # coverage = str_replace(coverage, "-Inf","0"),
         # coverage = as.numeric(coverage)
         ) %>%
  bed_to_granges ()



  mazia_musa_ac_align_cov_filtered_grange %>%
  tail


rm(mazia_musa_ac_align_cov)
rm (mazia_musa_sc_align_cov,mazia_ensete_gl_align_cov)

# Generate karyotype plot

library(karyoploteR)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")

library(GenomicRanges)

BiocManager::install("GenomicRanges")


# Genomic ranges 

MA_genome_bed %>% head()
# MB_genome_bed
# MS_genome_bed
# EG_genome_bed

## mazia vs musa ac

ymin = floor(min(mazia_musa_ac_align_cov_filtered_grange$id))
ymax = ceiling(max(mazia_musa_ac_align_cov_filtered_grange$id))


MA_genome_bed  %>%
  ggplot(aes(chr,end))+
  geom_point()
  # MA_genome_bed %>%
  # mutate(
  #   chr = str_remove(chr, "0hr"),
  #   chr = str_replace(chr,"chr0", "chr")) 

#Open the device
  MA_genome_bed_chr1_3 <-
  MA_genome_bed %>%
  dplyr::filter(chr == "chr01" | chr=="chr02" | chr =="chr03") %>%
  mutate(chr = str_remove(chr,"0"),
         chr = str_replace(chr, "chr","Chr"))

tiff("/Users/u1866313/rstudio/bgimazia_musa/Ensete_Musa/manuscript/mazia_musa_ac_chr1_3_karyoplot.tiff", width = 600, height = 500)

# mazia_musa_ac_align_cov_filtered_grange %>% GRanges(seqnames = "Chr1")
kp <- plotKaryotype(genome=MA_genome_bed %>%
                      dplyr::filter(chr == "chr01" | chr=="chr02" | chr =="chr03") %>%
                      mutate(chr = str_remove(chr,"0"),
                             chr = str_replace(chr, "chr","Chr")), ideogram.plotter = NULL, plot.type=6 )
kpAddCytobandsAsLine(kp)
kpBars(kp, mazia_musa_ac_align_cov_filtered_grange_chr1_3,  
       y1=mazia_musa_ac_align_cov_filtered_grange_chr1_3$id, 
       ymax=max(mazia_musa_ac_align_cov_filtered_grange_chr1_3$id)/2.5, 
       col="#006400", r0=0.5, r1=1, border=NA)
kpBars(kp, bedadti_musa_ac_align_cov_grange,  
       y1=bedadti_musa_ac_align_cov_grange$id, 
       ymax=max(bedadti_musa_ac_align_cov_grange$id)/2.5, 
       col="#0000FF", r0=0.5, r1=0,border=NA)

#Close the device
dev.off()


# mazia vs musa balbisiana





mazia_musa_ac_align_cov_filtered_grange_chr1_3$y <- mazia_musa_ac_align_cov_filtered_grange_chr1_3$id

ymin = floor(min(mazia_musa_ac_align_cov_filtered_grange_chr1_3$y))
ymax = ceiling(max(mazia_musa_ac_align_cov_filtered_grange_chr1_3$y))

kp <- plotKaryotype(genome=MA_genome_bed %>%
                      filter(chr == "chr01" | chr=="chr02" | chr =="chr03") %>%
                      mutate(chr = str_remove(chr,"0"),
                             chr = str_replace(chr, "chr","Chr")), ideogram.plotter = NULL, plot.type=6 )
kpPlotRegions(kp, data=mazia_musa_ac_align_cov_filtered_grange_chr1_3, 
              avoid.overlapping = FALSE, ymin=ymin, ymax=ymax, col="#006400",data.panel = 1)
kpPlotRegions(kp, bedadti_musa_ac_align_cov_grange,  
       ymin=ymin, 
       ymax=ymax, 
       col="#0000FF", data.panel = 1)



mazia_musa_ac_align_cov_filtered_grange$y <- mazia_musa_ac_align_cov_filtered_grange$id

ymin = floor(min(mazia_musa_ac_align_cov_filtered_grange$y))
ymax = ceiling(max(mazia_musa_ac_align_cov_filtered_grange$y))

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL,
                    plot.type=6, plot.params = pp, cex= 1.5,chromosomes = "chr01")
kpPlotRegions(kp, data=mazia_musa_ac_align_cov_filtered_grange_chr1_3,
              avoid.overlapping = FALSE, col="#006400",ymin=ymin, ymax=ymax,data.panel = 1)
kpPlotRegions(kp, data=bedadti_musa_ac_align_cov_grange,
              avoid.overlapping = FALSE, col="#0000FF",ymin=ymin, ymax=ymax,data.panel = 2)
# kpPlotRegions(kp, data=more.regions, data.panel=2, col="#0000FF") # Bedadeti

kpAxis(kp, data.panel=1, tick.pos = c(0, 3, 6))
kpAxis(kp, data.panel=2)

kp <- plotKaryotype(plot.type = 2)
kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6, chromosomes=c("Chr1", "Chr2") )
kpBars(kp, mazia_musa_ac_align_cov_filtered_grange,
       data.panel = 1,
       y1=mazia_musa_ac_align_cov_filtered_grange$id, 
       ymax = max(mazia_musa_ac_align_cov_filtered_grange$id), 
       col="blue")
kpBars(kp, 
       data=mazia_musa_ac_align_cov_filtered_grange, 
       data.panel = 2, 
       y1=mazia_musa_ac_align_cov_filtered_grange$id, 
       ymax=max(mazia_musa_ac_align_cov_filtered_grange$id), 
       col="red")



# bedadeti vs musa ac
ymin = floor(min(bedadti_musa_ac_align_cov_grange$id))
ymax = ceiling(max(bedadti_musa_ac_align_cov_grange$id))


bedadti_musa_ac_align_cov_grange %>%
  head()

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6, chromosomes = "c1")
kpBars(kp, bedadti_musa_ac_align_cov_grange,  y1=bedadti_musa_ac_align_cov_grange$y, ymax=max(bedadti_musa_ac_align_cov_grange$y), 
       col="#006400", border=NA)
kpAxis(kp, data.panel=1, tick.pos = c(0, 3, 6))
kpAxis(kp, data.panel=2)

kpLines(kp, data=bedadti_musa_ac_align_cov_grange, ymin=ymin, ymax=ymax)
kpAxis(kp, ymin = ymin, ymax=ymax)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="#EE82EE")
#######



bedadti_musa_ac_align_cov %>% tail()

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=2, plot.params = pp, cex= 1.5)
kpAddCytobandsAsLine(kp)
# kpAddMainTitle(kp, "M. acuminata", cex=2)
kpPlotRegions(kp, data=mazia_musa_ac_align_cov_filtered_grange_chr1_3, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=bedadti_musa_ac_align_cov_grange, avoid.overlapping = FALSE, col="red", data.panel=2)
# 
# kp <- plotKaryotype(plot.type = 6, genome = MA_genome_bed, main = "Coverage")
# kp <- kpPlotBAMCoverage(kp, data=bam1)
# 
# bam1 %>% head()
# 
# kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=4, plot.params = pp, cex= 1.5,chromosomes = "chr01")
# kpPlotRegions(kp, data=bedadti_musa_ac_align_cov, avoid.overlapping = FALSE, col="deepskyblue",ymin=0, ymax=coverage)


regions <- createRandomRegions(nregions=10000, length.mean = 1e6, mask=NA, non.overlapping=FALSE)

regions@ranges@start %>% tail()
library(BiocManager)
# install("BSgenome.Hsapiens.UCSC.hg19")
# BiocManager::install("plyranges")
# library(plyranges)

# bedadti_musa_ac_align_cov_filtered <-
# bedadti_musa_ac_align_cov %>% 
#   dplyr::rename(y = coverage) %>%
#   filter(y >= 1000) %>%
#   mutate(y = case_when(y >= 1000 ~ 1000)) %>%
#   rbind(
#     bedadti_musa_ac_align_cov %>% 
#       dplyr::rename(y = coverage) %>%
#       filter(y < 1000)
#   ) %>%
#   filter(str_detect(chr,"^chr")) %>%
#   select(chr,start, end,y)


bedadti_musa_ac_align_cov_filtered <-
  bedadti_musa_ac_align_cov %>% 
  dplyr::rename(y = coverage) %>% 
    filter(str_detect(chr,"^chr")) %>%
  mutate(y = round(log10(y),1),
         chr = str_remove(chr,"hr0"),
         chr = str_remove(chr, "hr"),
         y = str_replace(y, "-Inf","0"),
         y = as.numeric(y)) 

# bedadti_musa_ac_align_cov_filtered %>% 
#   mutate(y = str_replace(y, "-Inf","0"),
#          y = as.numeric(y)) %>% str()

# bedadti_musa_ac_align_cov_grange <- 
  bedadti_musa_ac_align_cov_filtered %>% head()
  bed_to_granges ()


rm(bedadti_musa_ac_align_cov_filtered,bedadti_musa_ac_align_cov)

bedadti_musa_ac_align_cov_grange %>%
  head()

bedadti_musa_ac_align_cov_grange$y <- bedadti_musa_ac_align_cov_grange$id

# kp <- plotKaryotype(chromosomes="chr01", plot.type=4)
# kpPlotHorizon(kp, data=rand.data, ymin=0, ymax=id)
ymin = floor(min(bedadti_musa_ac_align_cov_grange$id))
ymax = ceiling(max(bedadti_musa_ac_align_cov_grange$id))

# kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, 
#                     plot.type=4, plot.params = pp, cex= 1.5,chromosomes = "chr01")
# kpPlotRegions(kp, data=bedadti_musa_ac_align_cov_grange, 
#               avoid.overlapping = FALSE, col="deepskyblue",ymin=ymin, ymax=ymax)

bedadti_musa_ac_align_cov_grange %>% head()
  dplyr::filter(y <= 1)

MA_genome_bed <-
MA_genome_bed %>%
  mutate(
         chr = str_replace(chr, "c", "Chr")) 

#Open the device
tiff("mazia_musa_ac_ch1_2.coverage.tiff", width = 800, height = 500)

  
kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6)
kpLines(kp, data=bedadti_musa_ac_align_cov_grange, ymin=ymin, ymax=ymax)
kpAxis(kp, ymin = ymin, ymax=ymax)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="#EE82EE")

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6, chromosomes = "Chr1")
kpLines(kp, data=bedadti_musa_ac_align_cov_grange, ymin=ymin, ymax=ymax)
kpAxis(kp, ymin = ymin, ymax=ymax)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="#006400") # Mazia
kpPlotRegions(kp, data=more.regions, data.panel=2, col="#0000FF") # Bedadeti


## add centromer location 
# 
# myregoins <- toGRanges("cytoband.txt")
# 
# read.delim("cytoband.txt") %>%
#   bed_to_granges ()
# 
# 
# 
# library(karyoploteR)
# 
# cn.regs <- sort(createRandomRegions(nregions = 100))
# 
# cyto <- getCytobands(genome="hg19")
# #remove the unneded cytobands
# cyto <- cyto[overlapsAny(cyto, cn.regs)]
# 
# cyto %>% tail()
# 
# 
# cust.genome <- toGRanges(data.frame(c("1A", "1B", "4A"), c(1, 1,1), c(8e8, 8e8, 8e8)))
# myregions <- toGRanges("cytoband.txt")   
# kp <- plotKaryotype(genome=cust.genome)
# kpPlotRegions(kp, myregions)
# 
# kpText(kp, data=myregions, y=0.5, labels = myregions$gieStain, data.panel = "ideogram", cex=0.8)
# 
# kpPlotRegions(kp, data=myregions)

## 
mazia_musa_ac_align_cov %>% tail()
mazia_musa_ba_align_cov
mazia_musa_sc_align_cov
mazia_ensete_gl_align_cov


mazia_musa_ac_align_cov_filtered <-
  mazia_musa_ac_align_cov %>% 
  dplyr::rename(y = coverage) %>% 
  filter(str_detect(chr,"^chr")) %>%
mutate(y = round(log10(y),1),
       y = str_replace(y, "-Inf","0"),
       y = as.numeric(y)) 
  # %>%
  # mutate(y = round(log10(y),1),
  #        chr = str_remove(chr,"hr0"),
  #        chr = str_remove(chr, "hr"),
  #        y = str_replace(y, "-Inf","0"),
  #        y = as.numeric(y)) 

ceiling(max(mazia_musa_ac_align_cov_filtered$y))

## Generate grante 
# mazia_musa_ac_align_cov_filtered_grange <-
mazia_musa_ac_align_cov_filtered %>%
  bed_to_granges ()
  

## updated 

ymin = floor(min(bedadti_musa_ac_align_cov_grange$id))
ymax = ceiling(max(bedadti_musa_ac_align_cov_grange$id))


bedadti_musa_ac_align_cov_grange %>%
  head()

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6, chromosomes = "c1")
kpBars(kp, bedadti_musa_ac_align_cov_grange,  y1=bedadti_musa_ac_align_cov_grange$y, ymax=max(bedadti_musa_ac_align_cov_grange$y), 
       col="#006400", border=NA)
kpAxis(kp, data.panel=1, tick.pos = c(0, 3, 6))
kpAxis(kp, data.panel=2)

kpLines(kp, data=bedadti_musa_ac_align_cov_grange, ymin=ymin, ymax=ymax)
kpAxis(kp, ymin = ymin, ymax=ymax)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="#EE82EE")

## mazia vs musa ac

mazia_musa_ac_align_cov_filtered_grange %>% head()

MA_genome_bed <-
MA_genome_bed %>%
  mutate(
    chr = str_replace(chr, "c", "chr0"),
    chr = str_replace(chr,"chr010", "chr10"),
    chr = str_replace(chr,"chr011", "chr11"))  


kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6, chromosomes = "chr01")
kpBars(kp, mazia_musa_ac_align_cov_filtered_grange,  y1=mazia_musa_ac_align_cov_filtered_grange$id, ymax=max(mazia_musa_ac_align_cov_filtered_grange$id), 
       col="#006400", border=NA)
kpAxis(kp, data.panel=1, tick.pos = c(0, 3, 6))
kpAxis(kp, data.panel=2)






