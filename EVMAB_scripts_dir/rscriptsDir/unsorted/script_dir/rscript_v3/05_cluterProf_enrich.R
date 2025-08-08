# BiocManager::install("clusterProfiler")
library(tidyverse)
library(clusterProfiler)
library(ggtext)

# Global genes 
EV_MAB_custered_genes.v5 %>% 
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(count = n())

GO_terms_global_genes.v5 %>% 
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(count = n())

interproscan_GO_terms%>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(count = n())

emapper_GO_terms_retrieved %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(count = n())

proteinfer_GO_terms_retrieved %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(clust,gene_id) %>%
  distinct() %>%
  group_by(clust) %>%
  summarise(count = n())


# GO_terms_global_genes.v5 %>%
#   filter(str_detect(GO_terms,"GO")) %>%
#   # distinct(clust)
#   select(clust,gene_id) %>%
#   distinct() %>%
#   group_by(clust) %>%
#   summarise(count = n())


# create an empty vector for saving output

clustProfiler_GOenrichment_out <- c()

# For EV bedadeti

# select target gene set
target_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>%
  filter(clust == "EV_bedadeti_specific") %>%
  select(gene_id,GO_terms) %>%
  distinct()
  # filter(str_detect(GO_terms,"GO"))  
 
# Extract global genes set
GO_terms_global_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>%
    select(GO_terms, gene_id ) %>%
    distinct()

# List targe gene set
target_genes_list <-  unique(target_genes$gene_id)

# run clusterProfiler
enriched_out <-
  clusterProfiler::enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

if (length(enriched_out) != 0) {
  # convert it to dataframe and tidy-up result
  enriched_out_df <- enriched_out@result %>%
    #separate ratios into 2 columns of data
    separate(BgRatio,
             into = c("size.term", "size.category"),
             sep = "/") %>%
    separate(
      GeneRatio,
      into = c("size.overlap.term", "size.overlap.category"),
      sep = "/"
    ) %>%
    #convert to numeric
    mutate_at(
      vars(
        "size.term",
        "size.category",
        "size.overlap.term",
        "size.overlap.category"
      ),
      as.numeric
    ) %>%
    #Calculate k/K
    mutate("k.K" = size.overlap.term / size.term) %>%
    mutate(#frac_cov = CUTOFF,
      category = "EV bedadeti specific")
  
  
  ## add data
  
  # clustProfiler_GOenrich_out <-
  #   enriched_out_df  %>% # add it to your list
  #   rbind(clustProfiler_GOenrich_out)
  
  clustProfiler_GOenrichment_out <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrichment_out)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

row.names(clustProfiler_GOenrichment_out) = NULL

## plot 
# bedadeti_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  dplyr::filter(category == "EV bedadeti specific") %>%
  # dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
  # dplyr::filter(size.overlap.term > 1) %>%
  # dplyr::filter(p.adjust <= 0.09 ) %>%
  # dplyr::filter(p.adjust <= 0.06 ) %>%
  # filter(str_detect(description,"bidirectional"))
  select(GO_terms,p.adjust,description,size.overlap.term) %>%
  arrange(desc(size.overlap.term))
  mutate(description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated")) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )


## EV mazia specific 

# select target gene set
target_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>%
  filter(clust == "EV_mazia_specific") %>%
  select(gene_id,GO_terms) %>%
  distinct()

# filter(str_detect(GO_terms,"GO"))  

# Extract global genes set
GO_terms_global_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(GO_terms, gene_id ) %>%
  distinct()

# List target gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  clusterProfiler::enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

if (length(enriched_out) != 0) {
  # convert it to dataframe and tidy-up result
  enriched_out_df <- enriched_out@result %>%
    #separate ratios into 2 columns of data
    separate(BgRatio,
             into = c("size.term", "size.category"),
             sep = "/") %>%
    separate(
      GeneRatio,
      into = c("size.overlap.term", "size.overlap.category"),
      sep = "/"
    ) %>%
    #convert to numeric
    mutate_at(
      vars(
        "size.term",
        "size.category",
        "size.overlap.term",
        "size.overlap.category"
      ),
      as.numeric
    ) %>%
    #Calculate k/K
    mutate("k.K" = size.overlap.term / size.term) %>%
    mutate(#frac_cov = CUTOFF,
      category = "EV mazia specific")

  ## add data

  clustProfiler_GOenrichment_out <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrichment_out)
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrichment_out) = NULL

# plot
# mazia_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  dplyr::filter(category == "EV mazia specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
  # dplyr::filter(size.overlap.term > 1) %>%
  # dplyr::filter(p.adjust <= 0.09 ) %>%
  # dplyr::filter(p.adjust <= 0.06 ) %>%
  # filter(str_detect(description,"bidirectional"))
  select(GO_terms,p.adjust,description,size.overlap.term) %>%
  arrange(desc(size.overlap.term)) %>%
  filter(GO_terms== "GO:0005975")
  mutate(description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated")) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )

## MA specific 

# select target gene set
target_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>% 
  filter(clust == "MA_specific") %>%
  select(gene_id,GO_terms) %>%
  distinct()

# filter(str_detect(GO_terms,"GO"))  

# Extract global genes set
GO_terms_global_genes <-
  GO_terms_global_genes.v5 %>%
  filter(str_detect(GO_terms,"GO")) %>%
  select(GO_terms, gene_id ) %>%
  distinct()

# List target gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  clusterProfiler::enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

if (length(enriched_out) != 0) {
  # convert it to dataframe and tidy-up result
  enriched_out_df <- enriched_out@result %>%
    #separate ratios into 2 columns of data
    separate(BgRatio,
             into = c("size.term", "size.category"),
             sep = "/") %>%
    separate(
      GeneRatio,
      into = c("size.overlap.term", "size.overlap.category"),
      sep = "/"
    ) %>%
    #convert to numeric
    mutate_at(
      vars(
        "size.term",
        "size.category",
        "size.overlap.term",
        "size.overlap.category"
      ),
      as.numeric
    ) %>%
    #Calculate k/K
    mutate("k.K" = size.overlap.term / size.term) %>%
    mutate(#frac_cov = CUTOFF,
      category = "MA specific")
  
  
  ## add data
  
  # clustProfiler_GOenrich_out <-
  #   enriched_out_df  %>% # add it to your list
  #   rbind(clustProfiler_GOenrich_out)
  
  clustProfiler_GOenrichment_out <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrichment_out)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrichment_out) = NULL


# plot
# MA_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  distinct() %>%
  dplyr::filter(category == "MA specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
  # dplyr::filter(p.adjust < 0.65 ) %>%
  distinct() %>%
  select(GO_terms,p.adjust,description,size.overlap.term) %>%
  arrange(desc(size.overlap.term))  %>%
  filter(GO_terms== "GO:0005975")

  mutate(description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated")) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )



## MB specific 

# select target gene set
  target_genes <-
    GO_terms_global_genes.v5 %>%
    filter(str_detect(GO_terms,"GO")) %>% 
    filter(clust == "MB_specific") %>%
    select(gene_id,GO_terms) %>%
    distinct()
  
  # filter(str_detect(GO_terms,"GO"))  
  
  # Extract global genes set
  GO_terms_global_genes <-
    GO_terms_global_genes.v5 %>%
    filter(str_detect(GO_terms,"GO")) %>%
    select(GO_terms, gene_id ) %>%
    distinct()
  
  # List target gene set
  target_genes_list <-  unique(target_genes$gene_id)
  
  enriched_out <-
    clusterProfiler::enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

if (length(enriched_out) != 0) {
  # convert it to dataframe and tidy-up result
  enriched_out_df <- enriched_out@result %>%
    #separate ratios into 2 columns of data
    separate(BgRatio,
             into = c("size.term", "size.category"),
             sep = "/") %>%
    separate(
      GeneRatio,
      into = c("size.overlap.term", "size.overlap.category"),
      sep = "/"
    ) %>%
    #convert to numeric
    mutate_at(
      vars(
        "size.term",
        "size.category",
        "size.overlap.term",
        "size.overlap.category"
      ),
      as.numeric
    ) %>%
    #Calculate k/K
    mutate("k.K" = size.overlap.term / size.term) %>%
    mutate(#frac_cov = CUTOFF,
      category = "MB specific")
  
  
  ## add data
  
  # clustProfiler_GOenrich_out <-
  #   enriched_out_df  %>% # add it to your list
  #   rbind(clustProfiler_GOenrich_out)
  
  clustProfiler_GOenrichment_out <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrichment_out)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrichment_out) = NULL


# plot
# MB_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  distinct() %>%
  dplyr::filter(category == "MB specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>% 
  # dplyr::filter(p.adjust < 0.05 ) %>%
  select(GO_terms,p.adjust,description,size.overlap.term) %>%
  arrange(desc(size.overlap.term)) %>%
  filter(GO_terms== "GO:0005975")

  mutate(description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated")) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top",
    legend.justification='left'
    
  )





# ggsave(paste0(DATA_DIR,"plot_out/EV_mazia_specific_enrich.pdf"), width=5, height=4)
# ggsave(paste0(DATA_DIR,"plot_out/EV_mazia_specific_enrich.01.pdf"), width=5, height=4)
ggsave(paste0(DATA_DIR,"plot_out/EV_mazia_specific_enrich.allMapped.syn.02.pdf"), width=5, height=4)

## Bedadeti
# bedadeti_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  # clustProfiler_GOenrich_out_allMapped %>%
  dplyr::filter(category == "EV bedadeti specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  mutate(description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(
           description,
           "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
           "isopentenyl diphosphate biosynthetic process..."
         ),
         description = str_replace(description,"abscisic acid-activated signaling pathway","abscisic acid-activated<br>signaling pathway"),
         description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated")) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top",
    legend.justification='left'
  )

# bedadeti_enrich_plot
# ggsave(paste0(DATA_DIR,"plot_out/EV_bedadeti_specific_enrich.01.pdf"), width=5, height=4)
# ggsave(paste0(DATA_DIR,"plot_out/EV_bedadeti_specific_enrich.00.pdf"), width=6, height=12)
# ggsave(paste0(DATA_DIR,"plot_out/EV_bedadeti_specific_enrich.allMapped.pdf"), width=5, height=4)


mazia_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  # clustProfiler_GOenrich_out_allMapped %>%
  dplyr::filter(category == "EV mazia specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>% 
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  # distinct(description)
  
  mutate(description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"carbohydrate metabolic process","Carbohydrate<br>metabolic process"),
         description = str_replace(description,"transmembrane transport","Transmembrane<br>transport"),
         description = str_replace(description,"proteolysis","Proteolysis"),
         description = str_replace(description,"translation","Translation")
         ) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    legend.justification=c(-3,0)
  )

mazia_enrich_plot

ggsave(paste0(DATA_DIR,"plot_out/EV_mazia_specific_enrich.allMapped.pdf"), width=5, height=4)


#MA specific 

MA_enrich_plot <-
  clustProfiler_GOenrichment_out %>%
  dplyr::filter(category == "MA specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(aspects == 'Biological process') %>%
    dplyr::filter(p.adjust < 0.05 ) %>%
  # distinct(description)
  
  mutate(description = str_replace(description,"intracellular protein transport","Intracellular<br>protein transporty"),
         description = str_replace(description,"defense response","Defense response"),
         description = str_replace(description,"recognition of pollen","Recognition of pollen"),
         description = str_replace(description,"regulation of transcription, DNA-templated","Regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"protein phosphorylation","Protein phosphorylation")
         ) %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    legend.justification=c(-3,0)
  )

MA_enrich_plot
# ggsave(paste0(DATA_DIR,"plot_out/MA_specific_enrich.01.pdf"), width=5, height=4)
# 
# ggsave(paste0(DATA_DIR,"plot_out/MA_specific_enrich.pdf"), width=6, height=12)
# ggsave(paste0(DATA_DIR,"plot_out/MA_specific_enrich.allMapped.pdf"), width=5, height=4)
ggsave(paste0(DATA_DIR,"plot_out/MA_specific_enrich.allMapped.syn.02.pdf"), width=5, height=4)
mazia_enrich_plot | MA_enrich_plot|MB_enrich_plot




#MB specific 

MB_enrich_plot <-
clustProfiler_GOenrichment_out %>%
  dplyr::filter(category == "MB specific") %>%
  dplyr::select(-geneID, -Description) %>%
  dplyr::rename(GO_terms = ID) %>%
  left_join(go_odb_aspects_description) %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  
  mutate(aspects = str_replace_all(
    aspects,
    c(
      "biological_process" = "Biological process",
      "molecular_function" = "Molecular function",
      "cellular_component" = "Cellular component"
    )
  ),
  "log10_pvalue" = log10(p.adjust)) %>%
  dplyr::filter(aspects == 'Biological process') %>%
  # distinct(description)
  mutate(description = str_replace(description,"transmembrane transport","Transmembrane transport"),
         description = str_replace(description,"defense response","Defense response"),
         description = str_replace(description,"recognition of pollen","Recognition of pollen"),
         description = str_replace(description,"regulation of transcription, DNA-templated","Regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"protein phosphorylation","Protein phosphorylation")
  ) %>%
  dplyr::filter(description != "NA") %>%
  ggplot(aes(
    x = reorder(description, p.adjust),
    y = k.K,
    fill = log10_pvalue
  )) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = Count),
    hjust = 1.1,
    vjust = 0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_x_discrete(position = "top") +
  
  theme_classic() +
  coord_flip() +
  labs(
    y = "Proportion of significant genes over<br>the total genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    # legend.justification='right'
    legend.justification=c(0,0)
    
  )

MB_enrich_plot
# ggsave(paste0(DATA_DIR,"plot_out/MB_specific_enrich.01.pdf"), width=5, height=4)
# ggsave(paste0(DATA_DIR,"plot_out/MB_specific_enrich.pdf"), width=6, height=12)
# ggsave(paste0(DATA_DIR,"plot_out/MB_specific_enrich.allMapped.pdf"), width=5, height=4)
ggsave(paste0(DATA_DIR,"plot_out/MB_specific_enrich.allMapped.syn.02.pdf"), width=5, height=4)


library(patchwork)

mazia_enrich_plot | MA_enrich_plot|MB_enrich_plot
# ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.04.pdf"), width=11, height=6)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.02.pdf"), width=11, height=4)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.02.png"), width=11, height=4)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.03.pdf"), width=11, height=4.5)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.03.png"), width=11, height=4.5)

