library(karyoploteR)
library(regioneR)
library(zoo)

set.seed(1234)

#Parameters
data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")

num.data.points <- 3000
num.big.regions.up <- 30
num.big.regions.down <- 30

num.mid.regions <- 6000

num.marks <- 90

#Create the random fake data  

#Big regions
big.regs.up <- joinRegions(createRandomRegions(nregions = num.big.regions.up, length.mean = 20000000, length.sd = 10000000, non.overlapping = TRUE, mask=NA), min.dist = 1)
big.regs.down <- joinRegions(createRandomRegions(nregions = num.big.regions.down, length.mean = 10000000, length.sd = 5000000, non.overlapping = TRUE, mask=big.regs.up), min.dist = 1)

#Data points
data.points <- createRandomRegions(nregions = num.data.points, length.mean = 1, length.sd = 0, non.overlapping = TRUE, mask=NA)
mcols(data.points) <- data.frame(y=rnorm(n = num.data.points, 0.5, sd = 0.1))
dp.colors <- sample(head(data.points.colors, 2), size = num.data.points, replace = TRUE)

#and move the data points with the big regions
data.points[overlapsAny(data.points, big.regs.up)]$y <- data.points[overlapsAny(data.points, big.regs.up)]$y + runif(n=numOverlaps(data.points, big.regs.up), min = 0.1, max=0.3)
data.points[overlapsAny(data.points, big.regs.down)]$y <- data.points[overlapsAny(data.points, big.regs.down)]$y - runif(n=numOverlaps(data.points, big.regs.down), min = 0.1, max=0.3)

#markers
marks <- createRandomRegions(nregions = num.marks, length.mean = 1, length.sd = 0)
mcols(marks) <- data.frame(labels=paste0("rs", floor(runif(num.marks, min = 10000, max=99999))))

#medium regions
mid.regs <- createRandomRegions(nregions = num.mid.regions, length.mean = 5000000, length.sd = 1000000, non.overlapping = FALSE)


kp <- plotKaryotype(plot.type = 2, chromosomes = c("chr1", "chr2", "chr3"))

### Data Panel 1 ###

#Big regions
kpRect(kp, data = big.regs.up, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.8)
kpRect(kp, data = big.regs.down, y0=0, y1=1, col="#DDFFDD", border=NA, r0=0, r1=0.8)

#Data points
kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
kpPoints(kp, data=data.points, pch=16, cex=0.5, col=dp.colors, r0=0, r1=0.8)

#Mean and sd of the data points.  
for(chr in seqlevels(kp$genome)) {
  chr.dp <- sort(keepSeqlevels(x = data.points, value = chr, pruning.mode = "coarse"))
  rmean <- rollmean(chr.dp$y, k = 6, align = "center")  
  rsd <- rollapply(data = chr.dp$y, FUN=sd, width=6)
  kpLines(kp, chr = chr, x=start(chr.dp)[3:(length(chr.dp)-3)], y=rmean, col=data.points.colors[3], r0=0, r1=0.8)
  kpPlotRibbon(kp, chr=chr, data=chr.dp[3:(length(chr.dp)-3)], y0=rmean-rsd, y1=rmean+rsd, r0=0, r1=0.8, col="#FF336633", border=NA)
}

#Markers
kpPlotMarkers(kp, data=marks, label.color = "#333333", r1=1.1, cex=0.5, label.margin = 5)

### Data Panel 2 ###

#medium regions and their coverage

kpPlotRegions(kp, data = mid.regs, r0 = 0.2, r1=1, border=NA, data.panel=2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0, col=data.points.colors[2], data.panel = 2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0.12, col=data.points.colors[1], data.panel = 2)

kpText(kp, chr=seqlevels(kp$genome), y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, col="#444444", label="30x", cex=0.8, pos=2)
kpAbline(kp, h=0.4, data.panel = 2, r0=0.2, r1=0, col=data.points.colors[3])




library(plyranges)
# load library 
library(tidyverse)

library(karyoploteR)
library(GenomicRanges)
setwd("D:/Sadik/EV_MAB_coverage//")


bigwig_MZ_MA.bw <-
bigwig_MZ_MA.bw %>%
  as.data.frame() %>%
  mutate(y = log10(score),
         y = case_when(score=="-Inf" ~ 0,
                           TRUE ~ score)) %>%
  select(-score) %>%
  as_granges()

# bigwig_MZ_MA.bw

bigwig_BD_MA.file <- "bedadeti_Musa_acuminata.mapped.sorted.makrdup.bw"

bigwig_BD_MA.bw <-
  read_bigwig(bigwig_BD_MA.file) %>%
  as.data.frame() %>%
    mutate(score = log10(score),
           score = case_when(score=="-Inf" ~ 0,
                             TRUE ~ score)) %>%
  
  as_granges()

bigwig_MZ_MA.file <- "mazia_Musa_acuminata.mapped.sorted.makrdup.bw"

bigwig_MZ_MA.bw <- read_bigwig(bigwig_MZ_MA.file)



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

bigwig_BD_MA.bw <- 
bigwig_BD_MA.bw  %>% 
  as.data.frame() %>%
  filter(str_detect(seqnames,"chr")) %>%
  as_granges()


centro_density <-
  
  read.delim("ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
  bed_to_granges ()


kp <- plotKaryotype(genome=MAB_genome_bed,  plot.type=2 )
# kp <- plotKaryotype(genome=MAB_genome_bed, ideogram.plotter = NULL, plot.type=2 )
kpPlotCoverage(kp, data=bigwig_MZ_MA.bw, r0=0.5, r1=1, col="#006400", data.panel = 1)
bigwig_MZ_MA.file

kp <- kpPlotBigWig(kp, data=bigwig_MZ_MA.file, r0=0, r1=1.2, lwd=2, border = "#006400")
kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)
kp <- kpPlotBigWig(kp, data=bigwig_BD_MA.file, r0=0.25, r1=2.5, border="#0000FF", lwd=2)
kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.25, r1=0.45, border="red", lwd=2)





kp <- plotKaryotype(plot.type = 2, chromosomes = c("chr1", "chr2", "chr3"))

### Data Panel 1 ###

#Big regions
kpRect(kp, data = big.regs.up, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.8)
kpRect(kp, data = bigwig_MZ_MA.bw, y0=0, y1=1, col="#DDFFDD", border=NA, r0=0, r1=0.8)

#Data points
kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5,)
kpPoints(kp, bigwig_MZ_MA.bw, pch=16, cex=0.2, col="#006400", data.panel = 1)



kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )


kp <- kpPlotBigWig(kp, data=bigwig_MZ_MA.file, lwd=2, border = "#006400", data.panel = 1)
kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)
#

########

library(plyranges)



bigwig_BD_MA.file <- "bedadeti_Musa_acuminata.mapped.sorted.makrdup.bw"

bigwig_BD_MA.bw <- read_bigwig(bigwig_BD_MA.file)

bigwig_MZ_MA.file <- "mazia_Musa_acuminata.mapped.sorted.makrdup.bw"

bigwig_MZ_MA.bw <- read_bigwig(bigwig_MZ_MA.file)



bigwig_panEV_MA_1k = "panEV_Musa_acuminata.coverage.bw"
bigwig_panEV_MB_1k = "panEV_Musa_balbisiana.coverage.bw"


bigwig_panEV_MB_1k.v1 <-
read_bigwig(bigwig_panEV_MB_1k) %>%
  as.data.frame() %>%
  filter(str_detect(seqnames,"Bchr")) %>%
  mutate(seqnames = str_remove(seqnames,"B")) %>%
  as_granges() #%>%
  
export(bigwig_panEV_MB_1k.v1,con = "panEV_Musa_balbisiana.coverage.v1.bw", format = "bigWig")

export.bw("panEV_Musa_balbisiana.coverage.v1.bw")
  mutate()

panEV_Musa_balbisiana.coverage.bw


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


# custome_cytoband = read.table("https://raw.githubusercontent.com/bernatgel/karyoploter_examples/master/Examples/Tutorial/CustomGenomes/mygenome.txt", sep = "\t")

# plotKaryotype(genome = custome_cytoband)



kp <- plotKaryotype(genome=MAB_genome_bed,  ideogram.plotter = NULL, plot.type=2 )


kp <- kpPlotBigWig(kp, data=bigwig_panEV_MA_1k, border="green", lwd=3, data.panel = 1, r=0, r1=3)
kpBars(kp, centro_density_MA,
       y1=centro_density_MA$id,
       ymax=max(centro_density_MA$id),
       col="#555555", r0=-0.3, r1=0.1,border=NA)
kpBars(kp, centro_density_MB,
       y1=centro_density_MB$id,
       ymax=max(centro_density_MB$id),
       col="red", r0=-0.3, r1=-0.7,border=NA)
#
kp <- kpPlotBigWig(kp, data=bigwig_panEV_MB_1k, border="blue", lwd=3, data.panel = 2, r=0, r1=3)



# if (.Platform$OS.type != "windows") {
  bigwig.file <- system.file("extdata", "BRCA.genes.hg19.bw", package = "karyoploteR")
  brca.genes.file <- system.file("extdata", "BRCA.genes.hg19.txt", package = "karyoploteR")
  brca.genes <- toGRanges(brca.genes.file)
  seqlevelsStyle(brca.genes) <- "UCSC"
  
  kp <- plotKaryotype(zoom = brca.genes[1])
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0, r1=0.2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r1=0.2, data.panel = 1)
  kp <- kpPlotBigWig(kp, data=bigwig.file,  data.panel = 2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.25, r1=0.45, border="red", lwd=2)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.5, r1=0.7, ymin=0, ymax=1000, border="gold", col=NA)
  kpAxis(kp, r0=0.5, r1=0.7, ymin=0, ymax=1000)
  kp <- kpPlotBigWig(kp, data=bigwig.file, r0=0.75, r1=0.95, ymin=0, ymax="visible.region", border="orchid", col=NA)
  kpAxis(kp, r0=0.75, r1=0.95, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
# }
  
  
  
  
  bigwig_panEV_MA_1k = "panEV_Musa_acuminata.1k.coverage.bw"
  bigwig_panEV_MB_1k = "panEV_Musa_balbisiana.1k.coverage.bw"
  
  
  bigwig_panEV_MA_1k.cv1 <-
    read_bigwig(bigwig_panEV_MA_1k) %>%
    as.data.frame() %>%
    filter(str_detect(seqnames,"chr")) %>%
    mutate(score = log10(score),
           score = case_when(score=="-Inf" ~ 0,
                             TRUE ~ score)) %>%
    as_granges() #%>%  
  
  bigwig_panEV_MB_1k.v1 <-
    read_bigwig(bigwig_panEV_MB_1k) %>%
    as.data.frame() %>%
    filter(str_detect(seqnames,"Bchr")) %>%
    mutate(seqnames = str_remove(seqnames,"B")) %>%
    mutate(score = log10(score),
           score = case_when(score=="-Inf" ~ 0,
                             TRUE ~ score)) %>%
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
