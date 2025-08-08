# set



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
BiocManager::install("GenomicRanges")

install.packages("tidyverse")
BiocManager::install("ggbio")
BiocManager::install("plyranges")
install.packages("zoo")
library(BiocManager)
install("BSgenome.Hsapiens.UCSC.hg19")
