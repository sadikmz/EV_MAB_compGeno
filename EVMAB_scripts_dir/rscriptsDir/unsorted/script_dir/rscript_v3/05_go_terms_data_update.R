# BiocManager::install("clusterProfiler")
library(tidyverse)
library(clusterProfiler)
library(ggtext)
library(rcartocolor)


# Set of identified genes specific to each species/genotype 
EV_MAB_custered_genes.v5 %>%
  head()

EV_MAB_custered_genes.v5.1 %>%
  head()

# Join with GO terms annotation from interproscan

EV_MAB_custered_genes.v5 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      distinct(gene_id, GO_terms)
  ) %>%
  head()

interproscan_GO_terms <-
  EV_MAB_custered_genes.v5.1 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      distinct(gene_id, GO_terms)
  ) %>%
  # mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
  filter(str_detect(GO_terms,"GO")) #%>%
# head()


interproscan_GO_terms.v1 <-
  interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
    # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  # filter(str_detect(clust,"EV_m")) %>%
  select(gene_id,description,clust) %>% 
  # distinct() 
  group_by(clust,description) %>%
  summarise(total_genes = n(),
            gene_list = paste(sort(unique(gene_id)),collapse = ", ")) %>%
  as.data.frame() %>%
  arrange(desc(total_genes)) %>%
  as.data.frame()

