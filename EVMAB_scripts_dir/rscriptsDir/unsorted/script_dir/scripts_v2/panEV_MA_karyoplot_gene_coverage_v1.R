# set working directory 
# setwd("D:/Sadik/EV_MAB_coverage/")


# load library 
library(tidyverse)
library(regioneR)
library(zoo)
library(karyoploteR)
library(GenomicRanges)
library(plyranges)


#  Generate karyotype plot 


# read data 
## PanEV vs M. acuminata 
# panEV_MA_genomecov_gene <-
#   read.delim("bigwigs/panEV_Musa_acuminata.genomecov.gene.bed", header = F)

# panEV_MA_genomecov <-
#   read.delim("bigwigs/panEV_Musa_acuminata.genomecov.bed", header = F)
# 
# panEV_MA_genomecov_grange <- 
#   panEV_MA_genomecov %>%
#   dplyr::rename(seqnames = V1,
#                 start = V2,
#                 end = V3,
#                 coverage = V4) %>%
#   select(seqnames, start, end, coverage) %>%
#   dplyr::filter(str_detect(seqnames,"chr")) %>%
#   mutate(
#     #seqnames = str_remove(seqnames,"B"),
#     coverage = log10(coverage),
#     coverage = case_when(coverage=="-Inf" ~ 0,
#                       TRUE ~ coverage)) %>%
#   as_granges()



# rm(panEV_MA_genomecov)
## PanEV vs M. balbisiana 

# panEV_MB_genomecov_gene <-
#   read.delim("bigwigs/panEV_Musa_balbisiana.genomecov.gene.v1.bed", header = F)

# panEV_MB_genomecov <-
#   read.delim("bigwigs/panEV_Musa_balbisiana.genomecov.v1.bed", header = F)
# 
# 
# panEV_MB_genomecov_grange <-
#   panEV_MB_genomecov %>%
#   dplyr::rename(seqnames = V1,
#                 start = V2,
#                 end = V3,
#                 coverage = V4) %>%
#   select(seqnames, start, end, coverage) %>% 
#   filter(str_detect(seqnames,"Bchr")) %>%
#   mutate(seqnames = str_remove(seqnames,"B"),
#          coverage = log10(coverage),
#          coverage = case_when(coverage=="-Inf" ~ 0,
#                               TRUE ~ coverage)
#   ) %>%
#   dplyr::filter(str_detect(seqnames,"chr")) %>%
#   as_granges()


# Updated input data
panEV_MA_frac_ovl_syn %>%
  head()
# panEV_MB_frac_ovl_syn
# 
# MA_gene_bed
# MB_gene_bed 
# 
# MB_specific_genes
# MA_gene_complement
# MB_gene_complement


EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = round(max(fraction_overlap),2)) %>%
  ungroup() %>%
  # select(seqid, start, end, gene_id,ref_genome) %>%
  # group_by(ref_genome) %>%
  # summarise(count = n())
  filter(ref_genome == "musa_acuminata") %>%
  select(seqid, start, end, fraction_overlap) %>%
  # filter(fraction_overlap < 0.25) %>%
  distinct(seqid,start,end) %>%
  nrow()

# grange coverage 
panEV_MA_genomecov_grange <-
  EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = round(max(fraction_overlap),2)) %>%
  ungroup() %>%
  # select(seqid, start, end, gene_id,ref_genome) %>%
  # group_by(ref_genome) %>%
  # summarise(count = n())
  filter(ref_genome == "musa_acuminata") %>%
  select(seqid, start, end, fraction_overlap) %>%
    # head()
  # panEV_MA_frac_ovl_syn %>%
    dplyr::rename(seqnames = seqid,
                  coverage = fraction_overlap) %>%
    select(seqnames, start, end, coverage) %>%
    as_granges()


# MA genes assigned chromosomes

MA_gene_bed <-
EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query=="panEV_reads") %>%
  filter(ref_genome == "musa_acuminata") %>%
  filter(str_detect(seqid,"chr")) %>%
  distinct(seqid,start,end) 

MA_gene <- 
  MA_gene_bed %>% 
  filter(str_detect(seqid,"chr")) %>%
  dplyr::rename(seqnames = seqid) %>%
  as_granges()

# Repeat regions
MA_repeats <-
  read.delim(paste0(DATA_DIR,"Musa_acuminata.gene.complement.bed"), header = F) %>%
  filter(str_detect(V1,"chr")) %>%
  dplyr::rename(seqnames = V1,
                start = V2,
                end = V3) %>%
  as_granges()

# MA specific genes 
MA_specific_genes_grange <-
EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>%
  left_join(
    EV_MAB_custered_genes.v5.1 %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          distinct(gene_id, GO_terms)
      ) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          select(-GO_terms)
      ) %>%
      distinct() %>%
      
      filter(clust == "MA_specific" 
             # clust == "MB_specific" 
             # clust == "EV_mazia_specific" 
             # clust == "EV_bedadeti_specific"
      ) %>% 
      distinct(seqid, start, end, gene_id) %>%
      # coverage here is to show the position on the genes in the plot
      mutate(coverage = 1)) %>%
  mutate(coverage = str_replace_na(coverage, "NA"),
         coverage = as.numeric(coverage)) %>%
  filter(coverage != "NA") %>%
  dplyr::rename(seqnames = seqid) %>%
  dplyr::distinct(seqnames,start,end,coverage) #%>%
  # summary()
  # group_by(seqnames) %>%
  # summarise(count = n()) %>%
  # summarise(sum = sum(count))
  # as.data.frame() %>%
    # head()
  # as_granges()


# MA_specific_genes_grange <-
# MA_specific_genes %>% 
#   group_by(seqid,start,end) %>%
#   summarise(frac_ovrerlap_syn = max(frac_ovrerlap_syn)) %>%
#   as.data.frame() %>%
#   filter(str_detect(seqid,"chr")) %>%
#   distinct() %>%
#   dplyr::rename(seqnames = seqid,
#                 coverage = frac_ovrerlap_syn) %>%
#   # head()
#   as_granges()
  

# MB_specific_genes_grange <-
#   MB_specific_genes %>% 
#   group_by(seqid,start,end) %>%
#   summarise(frac_ovrerlap_syn = max(frac_ovrerlap_syn)) %>%
#   as.data.frame() %>%
#   filter(str_detect(seqid,"chr")) %>%
#   mutate(seqid = str_remove(seqid,"B")) %>% 
#   filter(str_detect(seqid,"chr")) %>%
#   distinct() %>%
#   dplyr::rename(seqnames = seqid,
#                 coverage = frac_ovrerlap_syn) %>%
#   # head()
#   as_granges()



# MB_gene <- 
# MB_gene_bed %>%
#   filter(str_detect(seqid,"chr")) %>%
#   mutate(seqid = str_remove(seqid,"B")) %>% 
#   dplyr::rename(seqnames = seqid) %>%
#   as_granges()





# MB_repeats <- 
#   MB_gene_complement %>% 
#   filter(str_detect(V1,"chr")) %>%
#   mutate(V1 = str_remove(V1,"B")) %>% 
#   dplyr::rename(seqnames = V1,
#                 start = V2,
#                 end = V3) %>%
#   as_granges()
# panEV_MB_genomecov <-
#   read.delim("bigwigs/panEV_Musa_balbisiana.genomecov.v1.bed", header = F)
# 
# 
# panEV_MB_genomecov_gene_grange <-
#   panEV_MB_genomecov %>%
#   dplyr::rename(seqnames = V1,
#                 start = V2,
#                 end = V3,
#                 coverage = V4) %>%
#   select(seqnames, start, end, coverage) %>% 
#   filter(str_detect(seqnames,"Bchr")) %>%
#   mutate(seqnames = str_remove(seqnames,"B"),
#          coverage = log10(coverage),
#          coverage = case_when(coverage=="-Inf" ~ 0,
#                               TRUE ~ coverage)
#   ) %>%
#   dplyr::filter(str_detect(seqnames,"chr")) %>%
#   as_granges()

# chromosome labesl 
MA_genome_bed <-
  read.delim(paste0(DATA_DIR,"EV_MAB_karyoplot/","Musa_acuminata.genome.bed")) %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(seqnames = seq.name) %>%
  as_granges()

# filter(chr == "chr01")


 
# filter(chr == "chr01")


# Species specific genes 

#  Chromosome length adjusted to both MA and MB
# MAB_genome_bed <-
#   MA_genome_bed %>%
#   left_join(Mb_genome_bed %>%
#               mutate(chr= str_remove (chr,"B")) %>%
#               dplyr::rename(end_MB = "end")) %>%
#   mutate(end_MAB = case_when(end_MB > end ~ end_MB,
#                              TRUE ~ end)) %>%
#   select(-end,-end_MB) %>%
#   dplyr::rename(end = end_MAB) %>%
#   bed_to_granges ()

# checking chromosome names 
# MAB_genome_bed <-
#   MAB_genome_bed %>%
#   as.data.frame() %>%
#   dplyr::filter(str_detect(seqnames,"chr")) %>%
#   as_granges()



# Centromere and karyotype range


# centro_density <-
#   read.delim(paste0(DATA_DIR,"EV_MAB_karyoplot/", "ALL_centro.density.100kb"), header= F, sep = ' ') %>%
#   dplyr::rename(chr = V1,
#                 start = V2,
#                 end = V3,
#                 density = V4) %>%
#   dplyr::filter(str_detect(chr,"^A")) %>%
#   mutate(chr = str_replace(chr, "A","chr")) %>%
#   # filter(chr == "chr01") %>%
#   bed_to_granges ()



centro_density_MA <-
  read.delim(paste0(DATA_DIR,"EV_MAB_karyoplot/", "ALL_centro.density.100kb"), header= F, sep = ' ') %>%
  dplyr::rename(chr = V1,
                start = V2,
                end = V3,
                density = V4) %>%
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>%
  bed_to_granges ()



# MB specific 

MB_specific_genes_grange <-
  EV_MAB_reads_coverage_geneID_GO_terms %>%
  filter(query == "panEV_reads") %>%
  filter(str_detect(ref_genome,"musa")) %>%
  distinct(seqid, start, end, gene_id,ref_genome, fraction_overlap) %>%
  group_by(seqid, start, end, gene_id, ref_genome) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>%
  left_join(
    EV_MAB_custered_genes.v5.1 %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          distinct(gene_id, GO_terms)
      ) %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms %>%
          select(-GO_terms)
      ) %>%
      distinct() %>%
      
      filter(clust == "MB_specific" 
             # clust == "MB_specific" 
             # clust == "EV_mazia_specific" 
             # clust == "EV_bedadeti_specific"
      ) %>% 
      distinct(seqid, start, end, gene_id) %>%
      # coverage here is to show the position on the genes in the plot
      mutate(coverage = 1)) %>%
  mutate(coverage = str_replace_na(coverage, "NA"),
         coverage = as.numeric(coverage)) %>%
  filter(coverage != "NA") %>%
  dplyr::rename(seqnames = seqid) %>%
  dplyr::distinct(seqnames,start,end,coverage) %>%
  group_by(seqnames) %>%
  # summarise(count = n()) %>%
  as.data.frame() %>%
  filter(str_detect(seqnames, "Bchr")) %>%
    mutate(seqnames = str_remove(seqnames, "B")) %>%
  as_granges()



# chr names

MB_genome_bed <-
  read.delim(paste0(DATA_DIR,"EV_MAB_karyoplot/","Musa_balbisiana.genome.bed")) %>%
  filter(str_detect(seq.name,"chr")) %>%
  dplyr::rename(seqnames = seq.name) %>%
  mutate(seqnames = str_remove(seqnames, "B")) %>%
  as_granges()
  # head()

centro_density_MB <- 
  read.delim(paste0(DATA_DIR,"EV_MAB_karyoplot/", "ALL_centro.density.100kb"), header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^B")) %>%
  mutate(chr = str_replace(chr, "B","chr")) %>% 
  bed_to_granges () 


ymax_uppper= ceiling(max(abs(range(panEV_MA_genomecov_grange$coverage))))
ymin_uppper= 1

ymax_lower= ceiling(max(abs(range(panEV_MB_genomecov_grange$coverage))))
ymin_lower= 1


## Open device 

# png(paste0(DATA_DIR,"EV_MAB_karyoplot/", "panEV_MAB_coverage_genomeWide.09.png"),width = 600, height = 500, units = "px")
png(paste0(result_dir, "panEV_MA_cov.png"),width = 800, height = 500, units = "px")

kp <- plotKaryotype(genome=MAB_genome_bed,  ideogram.plotter = NULL, plot.type=2, cex=1 )
# kpRect(kp, data = MA_gene, y0=0, y1=1, col="#5809eb", border=NA, r0=0, r1=0.2)
kpRect(kp, data = MA_gene, y0=0, y1=1, col="#5809eb", border=NA, r0=0, r1=0.3)
kpRect(kp, data = MA_repeats, y0=0, y1=1, col="#09ebd4", border=NA, r0=0, r1=-0.3)

# PanEV reads coverage against MA
kp <- kpBars(kp, MA_specific_genes_grange, 
             y1=MA_specific_genes_grange$coverage, 
               ymax = max(MA_specific_genes_grange$coverage),
             col="#8c8ca1", r0=0.31, r1=1.2,  border=NA )

kpAxis(kp, ymin = 0.25, ymax = 1, r0=0.21, r1=1.2, numticks = 3, col="#666666", cex=0.5)

# kpAxis(kp, ymax = ymax_uppper, ymin = ymin_uppper, r1 = 0.8, cex = 0.5 )

# MA centromeric repeats 

kpBars(kp, centro_density_MA,
       y1=centro_density_MA$id,
       ymax=max(centro_density_MA$id)/1.5,
       col="black", r0=-0.31, r1=-0.1,border=NA)

# MA specific genes 

# kp <- kpBars(kp, MA_specific_genes_grange, 
#              y1=MA_specific_genes_grange$coverage, 
#              ymax = max(MA_specific_genes_grange$coverage),
#              col="#fd0101", r0=-0.32, r1=-1.3, border=NA )

kp <- kpBars(kp, MB_specific_genes_grange, 
             y1=MB_specific_genes_grange$coverage, 
             ymax = max(MB_specific_genes_grange$coverage),
             col="#fd0101", r0=-0.32, r1=-1.3, border=NA )

# MB centromeric repeats
# kpBars(kp, centro_density_MB,
#        y1=centro_density_MB$id,
#        ymax=max(centro_density_MB$id)/1.5,
#        col="black", r0=-0.3, r1=-0.5,border=NA)

# 
# kp <- kpBars(kp, MA_specific_genes_grange, 
#              y1=MA_specific_genes_grange$coverage, 
#              ymax = max(MA_specific_genes_grange$coverage),
#              col="#0000FF", r0=-0.48, r1=-0.9, border=NA )

dev.off()

# kp <- kpBars(kp, panEV_MB_genomecov_grange, 
#              y1=panEV_MB_genomecov_grange$coverage, 
#              ymax = max(panEV_MB_genomecov_grange$coverage)/4,
#              col="#0000FF", r0=-0.48, r1=-0.9, border=NA )
# kpAxis(kp, ymax = ymax_uppper, ymin = ymin_uppper, r0  = 0.8, cex = 1 )


dev.off()



# centro_density_MA %>%
#   as.data.frame() %>%
#   filter(seqnames == "chr10") %>%
#   filter(start > 4700000) %>%
#   head()

# save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.05.RData"))



