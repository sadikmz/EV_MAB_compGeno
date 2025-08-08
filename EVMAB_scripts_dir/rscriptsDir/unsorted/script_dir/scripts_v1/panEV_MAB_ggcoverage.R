# load metadata
library(ggcoverage)
detach(package:ggcoverage, unload = TRUE)
# install via CRAN

# install via Github
# install.package("remotes")   #In case you have not installed it.
install.packages("remotes")
BiocManager::install("areyesq89/GenomeMatrix")
remotes::install_github("showteeth/ggcoverage")

library("rtracklayer")
library("ggcoverage")
library("ggpattern")

meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
sample.meta = read.csv(meta.file)
sample.meta
#>        SampleName    Type Group
#> 1 ERR127302_chr14 KO_rep1    KO
#> 2 ERR127303_chr14 KO_rep2    KO
#> 3 ERR127306_chr14 WT_rep1    WT
#> 4 ERR127307_chr14 WT_rep2    WT
#> 
#> 

# track folder
track.folder = system.file("extdata", "RNA-seq", package = "ggcoverage")
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         region = "chr14:21,677,306-21,737,601", extend = 2000,
                         meta.info = sample.meta)
# check data
head(track.df)
#>   seqnames    start      end score    Type Group
#> 1    chr14 21675306 21675950     0 KO_rep1    KO
#> 2    chr14 21675951 21676000     1 KO_rep1    KO
#> 3    chr14 21676001 21676100     2 KO_rep1    KO
#> 4    chr14 21676101 21676150     1 KO_rep1    KO
#> 5    chr14 21676151 21677100     0 KO_rep1    KO
#> 6    chr14 21677101 21677200     2 KO_rep1    KO
#> 
#> 
#> 

# create mark region
mark.region=data.frame(start=c(21678900,21732001,21737590),
                       end=c(21679900,21732400,21737650),
                       label=c("M1", "M2", "M3"))
# check data
mark.region
#>      start      end label
#> 1 21678900 21679900    M1
#> 2 21732001 21732400    M2
#> 3 21737590 21737650    M3

# Load gtf

gtf.file = system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf.gr = rtracklayer::import.gff(con = gtf.file, format = 'gtf')


# Basic coverage plot 



# sample.meta <- data.frame(
#   SampleName = c("Chr18_MCF7_ER_1", "Chr18_MCF7_ER_2", "Chr18_MCF7_ER_3", "Chr18_MCF7_input"),
#   Type = c("MCF7_ER_1", "MCF7_ER_2", "MCF7_ER_3", "MCF7_input"),
#   Group = c("IP", "IP", "IP", "Input")
# )
# # track folder
# track.folder <- system.file("extdata", "ChIP-seq", package = "ggcoverage")
# # load bigwig file
# track.df <- LoadTrackFile(
#   track.folder = track.folder, format = "bw",
#   meta.info = sample.meta
# )



# 
# joint view
# Create line plot for every sample (facet.key = "Type") and color by every sample (group.key = "Type"):
  
basic.coverage = ggcoverage(data = track.df, color = "auto", 
                              plot.type = "joint", facet.key = "Type", group.key = "Type",
                              mark.region = mark.region, range.position = "out")
basic.coverage



# Create group average line plot (sample is indicated by facet.key = "Type", group is indicated by group.key = "Group"):

basic.coverage = ggcoverage(data = track.df, 
                            color = "auto", 
                            plot.type = "joint", 
                            facet.key = "Type", 
                            group.key = "Group", 
                            joint.avg = TRUE,
                            mark.region = mark.region, 
                            range.position = "out")
basic.coverage



# facet view
basic.coverage = ggcoverage(data = track.df, 
                            color = "auto", 
                            plot.type = "facet",
                            mark.region = mark.region, 
                            range.position = "out")
basic.coverage







# Custom Y-axis style
# Change the Y-axis scale label in/out of plot region with range.position:
  
library(tidyverse)
library(plyranges)
track.df %>%
head()
    basic.coverage = ggcoverage(data = track.df, 
                              color = "auto", 
                              plot.type = "facet",
                              mark.region = mark.region, 
                              group.key = "Group",
                              range.position = "in")
basic.coverage



########


track.folder = "D:/Sadik/EV_MAB_coverage/test//"
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         # region = "chr14:21,677,306-21,737,601", extend = 2000,
                         meta.info = sample.meta)




read_bigwig("C:/Users/Sadik.DESKTOP-UQ756TB/AppData/Local/R/win-library/4.3/ggcoverage/extdata/RNA-seq/ERR127302_chr14.bw")

panEV_MA_cov_bw_grange <-
read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Arkiya_362.bw") %>%
  as.data.frame() %>% 
  mutate(Type = "Arkiya_362",
         Group = "Domesticated\nenset") %>%

bind_rows(
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Arkiya_362.bw") %>%
      as.data.frame() %>% 
      mutate(Type = "Arkiya_455",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Astara.bw") %>%
      as.data.frame() %>% 
      mutate(Type = "Astara",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Bedadeti.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Bedadeti",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Buffero.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Buffero",
             Group = "Domesticated\nenset"),
    #   
      read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Derea.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Derea",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Erpha13.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Erpha13",
             Group = "Wild\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Erpha20.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Erpha20",
             Group = "Wild\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Lochinge_253.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Lochinge_253",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Mazia_208.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Mazia_208",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Nechuwe.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Nechuwe",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Nobo.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Nobo",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Onjamo.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Onjamo",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Siyute.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Siyute",
             Group = "Domesticated\nenset"),

    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Yako.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Yako",
             Group = "Domesticated\nenset"),


    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.ma/Yanbule.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Yanbule",
             Group = "Domesticated\nenset")
    
  )
  

  

MA_gene_gff3 =rtracklayer::import.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musa_acuminata_pahang_v4.gff3",
                                     format = "gff3")

rtracklayer::export(MA_gene_gff3, "D:/Sadik/EV_MAB_coverage/coverage_noreads/Musa_acuminata_pahang_v4.gtf")

MA_gene_gtf = rtracklayer::import.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musa_acuminata_pahang_v4.gtf")


library(ape)
MA_gene_gtf = read_gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musa_acuminata_pahang_v4.gff3")

# MA_gene_gtf<- ape::read.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musa_acuminata_pahang_v4.gff3") #%>%
  # head()

MA_gene_gtf %>%
head()




panEV_MA_cov_bw_grange_chr10 <-
  panEV_MA_cov_bw_grange %>%
  filter(seqnames == "chr10") %>%
  # filter(Type == "Arkiya_362" ) %>%
  mutate(#seqnames = str_remove(seqnames,"B"),
    score = log10(score),
    score = case_when(score=="-Inf" ~ 0,
                      TRUE ~ score)) %>%
  as.data.frame() %>%
  filter(start >= 4322145) %>%
  filter(end <= 6230000) 
# head()


# set chromosome and mark region

mark_region <-
MA_gene_gtf %>%
  filter(seqnames == "chr10") %>%
  filter(type =="gene") %>%
  as.data.frame() %>%
  filter(start >= 4322145) %>%
  filter(end <= 6230000) #%>%
  # head()
  # as_granges()
  # select(seqnames,start,end) 


# # mark_region <-
#   MA_gene_gtf %>%
#   filter(seqnames == "chr10") %>%
#   # filter(type =="gene") %>%
#   as.data.frame() %>%
#   filter(start >= 4322145) %>%
#   filter(end <= 6230000) %>%
#     mutate(Note = str_replace_na(Note,"NA")) %>%
#     nrow()
#     distinct(Note)
# # head()



mark_region_gtf <-
  MA_gene_gtf %>%
  filter(seqnames == "chr10") %>%
  # filter(type =="gene") %>%
  # as.data.frame() %>%
  filter(start >= 4322145) %>%
  filter(end <= 6230000) %>%
    # dplyr::rename(seqnames = seqid) %>%
  # head()
    # as.data.frame() %>% 
    mutate(gene_type = "UNKKOWN",
           gene_name = Parent ) #%>%
    # as_granges()

# 4,322,154-6,234,365


# plot

# facet view
basic.coverage = ggcoverage(data = panEV_MA_cov_bw_grange_chr10,
                            # color = "auto", 
                            plot.type = "facet",
                            mark.region = mark_region,
                            range.position = "out",
                            # range.size =200000,
                            # mark.label.size=10,
                            joint.avg = TRUE,
                            mark.alpha =0.01,
                            facet.key = "Group")

basic.coverage

basic.coverage +
  geom_gene(gtf.gr = mark_region_gtf,
            plot.height = 0.2,
            overlap.style = "tight",
            plot.space = 0.01,
            # arrow.num = 10,
            gene.size = 4,
            exon.size=1,
            arrow.size = 0
            )

ggsave("D:/Sadik/EV_MAB_coverage/EV_MA_syntenic_region.png", width = 5, height = 3)


meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
sample.meta <- utils::read.csv(meta.file)
# track folder
track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
# load bigwig file
track.df <- LoadTrackFile(track.folder = track.folder, format = "bw",region = "chr14:21,677,306-21,737,601",
                          extend = 2000, meta.info = sample.meta)
gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
ggcoverage(data = track.df, color = "auto", range.position = "out")


# MA specific regions 
# chr11:13,429,867-14,630,589
chr07:779,201-797,326
chr07:38,308,582-38,329,630
chr07:38,059,796-38,091,476
chr07:36,858,009-36,877,578
# 4,322,154-6,234,365
# poorly aligned 
chr02:1-15,465,198
chr02:158,690-243,227
168,253-174,997
chr02:1-412,852
# 4,322,154-6,234,365
chr02:2,133,199-2,147,796
chr02:2,151,597-2,254,209
chr02:2,152,147-2,290,464
chr02:2,684,614-2,801,218
chr02:3,403,925-4,231,821
chr01:24,372,303-30,726,606


# Telomeric Musa specific 
chr06:1-75,972
# 43,117,341
chr06:43,094,134-43,098,429
chr05:1-176,111
chr05:428,740-473,091
chr05:3,306,416-3,414,751
chr04:42,528,799-42,648,559
chr04:42,496,129-42,781,943
chr03:40,532,431-40,587,160#
chr03:41,185,616-41,790,650

panEV_MA_cov_bw_grange_chr02 <-
  panEV_MA_cov_bw_grange %>%
  filter(seqnames == "chr03") %>%
  # filter(Type == "Arkiya_362" ) %>%
  mutate(#seqnames = str_remove(seqnames,"B"),
    score = log10(score),
    score = case_when(score=="-Inf" ~ 0,
                      TRUE ~ score)) %>%
  as.data.frame() %>%
  filter(start >= 40532431) %>%
  filter(end <= 40587160) 
# head()+1


mark_region <-
  MA_gene_gtf %>%
  filter(seqnames == "chr03") %>%
  filter(type =="gene") %>%
  as.data.frame() %>%
  filter(start >= 40532431) %>%
  filter(end <= 40587160)  %>%
      mutate(Note = str_replace_na(Note,"NA"))



mark_region_gtf <-
  MA_gene_gtf %>%
  filter(seqnames=="chr03") %>%
  # filter(type =="gene") %>%
  filter(start >= 40532431) %>%
  filter(end <= 40587160)  %>%
  mutate(gene_type = "UNKKOWN",
         gene_name = Parent ) #%>%


# plot

# facet view


basic.coverage = ggcoverage(data = panEV_MA_cov_bw_grange_chr02,
                            # color = "auto", 
                            plot.type = "facet",
                            mark.region = mark_region,
                            range.position = "out",
                            # range.size =200000,
                            # mark.label.size=10,
                            joint.avg = TRUE,
                            mark.alpha =0.01,
                            facet.key = "Group")

basic.coverage

basic.coverage +
  geom_gene(gtf.gr = mark_region_gtf,
            plot.height = 0.2,
            overlap.style = "tight",
            plot.space = 0.01,
            arrow.num = 20,
            gene.size = 1.5,
            exon.size=1,
            arrow.size = 1
  )

ggsave("D:/Sadik/EV_MAB_coverage/EV_MA_syntenic_region.PAV.telomeric.chr03.00.png", width = 5, height = 3)





##### MB

panEV_MB_cov_bw_grange <-
  read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Arkiya_362.bw") %>%
  as.data.frame() %>% 
  mutate(Type = "Arkiya_362",
         Group = "Domesticated\nenset") %>%
  
  bind_rows(
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Arkiya_362.bw") %>%
      as.data.frame() %>% 
      mutate(Type = "Arkiya_455",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Astara.bw") %>%
      as.data.frame() %>% 
      mutate(Type = "Astara",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Bedadeti.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Bedadeti",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Buffero.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Buffero",
             Group = "Domesticated\nenset"),
    #   
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Derea.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Derea",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Erpha13.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Erpha13",
             Group = "Wild\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Erpha20.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Erpha20",
             Group = "Wild\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Lochinge_253.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Lochinge_253",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Mazia_208.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Mazia_208",
             Group = "Domesticated\nenset"),
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Mazia_stlfr.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Mazia_stlfr",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Nechuwe.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Nechuwe",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Nobo.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Nobo",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Onjamo.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Onjamo",
             Group = "Domesticated\nenset"),
    # 
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Siyuti.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Siyute",
             Group = "Domesticated\nenset"),
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Yako.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Yako",
             Group = "Domesticated\nenset"),
    
    
    read_bigwig("D:/Sadik/EV_MAB_coverage/track.folder.gene.mb/Yanbule.bw") %>%
      as.data.frame() %>%
      mutate(Type = "Yanbule",
             Group = "Domesticated\nenset")
    
  )

MB_gene_gff3 =rtracklayer::import.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musba.gff3",
                                      format = "gff3")

rtracklayer::export(MB_gene_gff3, "D:/Sadik/EV_MAB_coverage/coverage_noreads/Musba.gtf")

MB_gene_gtf = rtracklayer::import.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musba.gtf")



# MB_gene_gtf = ape::read.gff("D:/Sadik/EV_MAB_coverage/coverage_noreads/Musba.gff3") %>%
  # dplyr::rename(seqnames = seqid) 
  

#centromeri region 
chr04:12,337,674-21,222,081
chr03:12,487,846-23,914,006
chr04:11,938,988-23,716,321
chr01:20,630,371-26,990,834




Bchr01:2,301,260-2,427,588
Bchr01:2,299,057-2,427,969
Bchr01:4,435,583-4,482,224#
Bchr01:19,507,961-19,559,906
Bchr01:21,084,701-21,269,638
Bchr01:33,709,878-35,342,771
Bchr02:23,896,732-23,976,268 #
Bchr10:12,207,926-13,128,327 # Centromeric 
Bchr11:313,066-348,872
Bchr11:387,327-432,050
Bchr11:3,717,514-3,743,361
Bchr10:39,847,308-39,886,822

panEV_MB_cov_bw_grange_chr02 <-
  panEV_MB_cov_bw_grange %>%
  filter(seqnames == "Bchr10") %>%
  mutate(#seqnames = str_remove(seqnames,"B"),
    score = log10(score),
    score = case_when(score=="-Inf" ~ 0,
                      TRUE ~ score)) %>%
  as.data.frame() %>%
  filter(start >= 39845000) %>%
  filter(end <= 39887000)  


# Mark region
mark_region <-
  MB_gene_gtf %>%
  filter(seqnames == "Bchr10") %>%
  filter(type =="gene") %>%
  as.data.frame() %>%
  filter(start >= 39845000) %>%
  filter(end <= 39887000)       # %>%
      # mutate(Note = str_replace_na(Note,""),
      #        label = Note)
  


# gene regions 
mark_region_gtf <-
  MB_gene_gtf %>%
  filter(seqnames=="Bchr10") %>%
  # filter(type =="gene") %>%
  filter(start >= 39845000) %>%
  filter(end <= 39887000)        %>%
  mutate(gene_type = "UNKKOWN",
         gene_name = Parent,
         Note = str_replace_na(Note,""),
         label = Note)



# plot

# facet view

basic.coverage = ggcoverage(data = panEV_MB_cov_bw_grange_chr02,
                            # color = "auto", 
                            plot.type = "facet",
                            mark.region = mark_region,
                            range.position = "out",
                            # range.size =200000,
                            # mark.label.size=1,
                            joint.avg = TRUE,
                            mark.alpha =0.01,
                            facet.key = "Group")

basic.coverage

basic.coverage +
  geom_gene(gtf.gr = mark_region_gtf,
            plot.height = 0.2,
            overlap.style = "tight",
            plot.space =5,
            arrow.num = 5,
            gene.size = 3,
            exon.size=1,
            arrow.size = 1,
            # overlap.gene.gap = 1
  )

# ggsave("D:/Sadik/EV_MAB_coverage/EV_MB.PAV.chr01.01.png", width = 5, height = 3)
# filter(start >= 4440000) %>%
#   filter(end <= 4483000) 
# ggsave("D:/Sadik/EV_MAB_coverage/EV_MB.PAV.chr01.02.png", width = 5, height = 3)
# filter(start >= 24000000) %>%
#   filter(end <= 25000000)  
# ggsave("D:/Sadik/EV_MAB_coverage/EV_MB.PAV.chr01.centromeric.png", width = 5, height = 3)
# filter(start >= 12700000) %>%
#   filter(end <= 14500000)


# ggsave("D:/Sadik/EV_MAB_coverage/EV_MB.PAV.chr10.00.png", width = 5, height = 3)
# filter(start >= 39845000) %>%
#   filter(end <= 39887000) 
save.image("D:/Sadik/EV_MAB_coverage/panEV_MB_coverage.RData")












