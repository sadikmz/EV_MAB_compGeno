# BiocManager::install("clusterProfiler")
library(tidyverse)
library(clusterProfiler)
library(ggtext)

# Set of identified genes specific to each species/genotype 

EV_MAB_custered_genes.v5.1 %>%
  head()


# emapper_GO_terms <-

emapper_annotations <-
read.delim(paste0(DATA_DIR, "eggnog_mapper_mazia_AED25.out/eggNOG.emapper.annotations"),
           skip = 4,
           sep = '\t') %>%
  rbind(
    read.delim(paste0(DATA_DIR,"eggnog_mapper_bedadeti_AED25.out/eggNOG.emapper.annotations"),
      skip = 4,
      sep = '\t'
    ),
    read.delim(paste0(DATA_DIR,"eggnog_mapper_musa_ac.out/eggNOG.emapper.annotations"),
      skip = 4,
      sep = '\t'
    ),
    read.delim(paste0(DATA_DIR,"eggnog_mapper_musa_ba.out/eggNOG.emapper.annotations"),
      skip = 4,
      sep = '\t'
    )
  ) 


# Generate a list of column names to seprate GO-term values

## End of column 

col_range <-
  
  emapper_annotations %>%
  rename(gene_id = X.query) %>%
  select(gene_id, GOs) %>%
  filter(str_detect(GOs, "GO:")) %>%
  mutate(count_sep = str_count(GOs, ",")) %>%
  arrange(desc(count_sep)) %>%
  # head()
  select(count_sep) %>%
max()

# col_range=135

col_range <- 
  seq(col_range) %>% 
  as_tibble() %>% 
  mutate(cols = str_c("A",value)) %>%
  select(cols) 

col_range = col_range$cols

start_pos = col_range %>% head()
end_pos = col_range %>% tail()

# extracted emapepr GO-terms 
emapper_GO_terms <-
emapper_annotations %>%
  rename(gene_id = X.query) %>%
  select(gene_id, GOs) %>%
  filter(str_detect(GOs, "GO:")) %>%
  mutate(count_sep = str_count(GOs, ",")) %>%
  arrange(desc(count_sep)) %>%
  # separate(GOs, sprintf("%s%02d", "A", col_range), sep = "\\,") %>%
  separate(GOs, into=col_range, sep = "\\,") %>%
  # head() %>%
  select(gene_id, start_pos[1]:end_pos[6]) %>%
  pivot_longer(!gene_id, names_to = "names", values_to = "GO_terms") %>%
  as.data.frame() %>%
  select(-names) %>%
  distinct

# GO_terms retrived from proteinfer annotation 
# EV_MAB_protinfer <-
#   read.delim(paste0("~/Downloads/", "mazia.protinfer.predictions.tsv"), 
#              header = T) %>%
#   rbind(
#     read.delim(paste0("~/Downloads/",  "bedadeti.protinfer.predictions.tsv"), header = T),
#     read.delim(paste0("~/Downloads/",  "musa_ac.protinfer.predictions.tsv"), header = T),
#     read.delim(paste0("~/Downloads/",  "musa_ba.protinfer.predictions.tsv"), header = T))

emapper_GO_terms_retrieved <-
  EV_MAB_custered_genes.v5.1 %>%
  left_join(EV_MAB_reads_coverage_geneID_GO_terms %>%
              distinct(gene_id, GO_terms)) %>%
  filter(!str_detect(GO_terms, "GO")) %>%
  distinct() %>%
  mutate(gene_list = "nogo") %>%
  select(-GO_terms) %>%
  # # find matching genes from no-go list
  left_join(
    emapper_GO_terms %>%
      mutate(
        gene_id = str_remove(gene_id, "-R\\w"),
        GO_terms = str_replace_na(GO_terms, "NA")
      ) %>%
      filter(GO_terms != "NA") %>%
      distinct()
  ) %>%
  mutate(GO_terms = str_replace_na(GO_terms, "NA")) %>%
  distinct()

emapper_GO_terms_retrieved %>%
  filter(gene_id == "Macma4_01_g01750.1") %>%
  head()

# Proteinfer reterived GO-terms
proteinfer_GO_terms_retrieved <-
  emapper_GO_terms_retrieved %>%
  filter(GO_terms == 'NA') %>%
  select(-GO_terms) %>%
  left_join(
    EV_MAB_protinfer %>%
      rename(gene_id = sequence_name) %>%
      filter(str_detect(predicted_label, "GO:")) %>%
      rename(GO_terms = predicted_label) %>%
      left_join(go_odb_aspects_description) %>%
      # filter(seq.name =="EVMZ.1.000227-RA") %>%
      select(-defination) %>%
      filter(aspects == "biological_process") %>%
      # filter(!str_detect(seq.name,"EVMZ")) %>%
      filter(confidence > 0.9) %>%
      separate(
        gene_id,
        into = c("gene_id", "isoform"),
        sep = "\\-"
      ) %>%
      select(gene_id, GO_terms) %>%
      distinct()
  ) %>%
  mutate (GO_terms = str_replace_na(GO_terms, "NA")) %>%
  # filter(GO_terms != "NA") %>%
  distinct()


proteinfer_GO_terms_retrieved %>%
  filter(gene_id == "Macma4_01_g01750.1") %>%
  head()


# Global genes 

GO_terms_global_genes.v5 <-
  interproscan_GO_terms #%>%
  # distinct() %>%
  # rbind(
  #   emapper_GO_terms_retrieved %>%
  #     select(-gene_list) %>%
  #     distinct(),
  #   proteinfer_GO_terms_retrieved %>%
  #     select(-gene_list) %>%
  #     distinct()
  #   )

  
GO_terms_global_genes.v5 %>%
    distinct() %>%
  # rbind(
  #     emapper_GO_terms_retrieved %>%
  #       select(-gene_list) %>%
  #       distinct()#,
  #     # proteinfer_GO_terms_retrieved %>%
  #     #   select(-gene_list) %>%
  #     #   distinct()
  #   ) %>%
  # filter(gene_id == "Macma4_01_g01750.1") %>%
    filter(str_detect(GO_terms,"GO")) %>%
    select(clust,gene_id) %>%
    # filter(GO_terms == "GO:0005975")
    distinct() %>%
    group_by(clust) %>%
    summarise(count = n())


EV_MAB_custered_genes.v5 %>% 
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(total = n()) %>%
  left_join(
interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(interproscan = n())) %>%
  left_join(

emapper_GO_terms_retrieved %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(emapper = n())) %>%
  left_join(

proteinfer_GO_terms_retrieved %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(proteinfer = n())) %>%
  mutate(
    total_with_goterms = interproscan + emapper +proteinfer,
    per = round(total_with_goterms/total * 100,2),
    tot = str_c(total_with_goterms," ", "(", per, "%",")")
  ) %>%
  write.table(paste0(result_dir,"identified_genes_with_go_terms_summary.txt"),
              col.names = T,
              row.names = F,
              quote = F,
              sep = "\t")


