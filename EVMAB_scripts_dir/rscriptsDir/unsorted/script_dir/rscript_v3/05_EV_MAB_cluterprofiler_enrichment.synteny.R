# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(ggtext)

# EV-MAB gene set
# # EV bedadeti gene sets
# EV_bedadeti_specific
# EV_bedadeti_MAB_core
# EV_bedadeti_MA_specific
# EV_bedadeti_MB_specific
# EV_bedadeti_noEVReadSupport_Butsyntenic2MAB
# 
# # EV mazia gene sets
# EV_mazia_specific
# EV_mazia_MAB_core
# EV_mazia_MA_specific
# EV_mazia_MB_specific
# EV_mazia_noEVReadSupport_Butsyntenic2MAB
# 
# # Musa acuminata gene sets
# MA_specific
# MA_EV_core
# MA_noMAReadSupport_Butsyntenic2EV
# 
# # Musa balbisiana gene sets
# MB_specific
# MB_EV_core
# MB_noMBReadSupport_Butsyntenic2EV

## update/add GO_terms annotation from genes who does not have InterProScan GO_terms annotation using eggNOG and proteinfer annotation
## Genes without interproscan GO_terms annotation 
EV_mazia_specific %>%
  rbind(
    EV_bedadeti_specific,
    MA_specific,
    MB_specific
  ) %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  distinct() %>%
  head()

max_col_range=163

# GO_terms retrived from eggNOG annotation 
eggnog_MAB_retrieved <- 
read.delim("/Users/u1866313/Downloads/out.emapper.annotations",skip = 4, sep = '\t') %>%
  rename(gene_id = X.query) %>%
  select(gene_id,GOs) %>%
  filter(str_detect(GOs,"GO:")) %>%
  mutate(count_sep = str_count(GOs,",")) %>%
  arrange(desc(count_sep)) %>%
  separate(GOs,sprintf("%s%02d", "A", col_range), sep = "\\,") %>%
  select(gene_id,A01:A163) %>%
  # head()
  # pivot_longer( cols = c("A01":"A163"), names_to = "names", values_to = "GO_terms") %>%
  pivot_longer( !gene_id, names_to = "names", values_to = "GO_terms") %>%
  as.data.frame() %>%
  select(-names) %>%
  distinct


# GO_terms retrived from proteinfer annotation 
EV_MAB_protinfer <-
  read.delim(paste0("~/Downloads/", "mazia.protinfer.predictions.tsv"), 
             header = T) %>%
  rbind(
    read.delim(paste0("~/Downloads/",  "bedadeti.protinfer.predictions.tsv"), header = T),
    read.delim(paste0("~/Downloads/",  "musa_ac.protinfer.predictions.tsv"), header = T),
    read.delim(paste0("~/Downloads/",  "musa_ba.protinfer.predictions.tsv"), header = T))

EV_MAB_protinfer_retrieved <- 
EV_mazia_specific %>%
  rbind(
    EV_bedadeti_specific,
    MA_specific,
    MB_specific
  ) %>%
  filter(!str_detect(GO_terms,"GO")) %>%
  distinct() %>%
  mutate(gene_list="nogo") %>%
  select(-GO_terms) %>%
  # filter(str_detect(gene_id,"EV")) %>%
  # head()
  # find matching genes from no-go list
  left_join(
    EV_MAB_protinfer %>%
      rename(gene_id = sequence_name) %>%
      filter(str_detect(predicted_label,"GO:")) %>%
      rename(GO_terms = predicted_label) %>%
      left_join(go_odb_aspects_description) %>%
      # filter(seq.name =="EVMZ.1.000227-RA") %>%
      select(-defination) %>%
      filter(aspects == "biological_process") %>%
      # filter(!str_detect(seq.name,"EVMZ")) %>%
      filter(confidence > 0.9 ) %>%
      separate(gene_id,into = c("gene_id","isoform"),sep = "\\-") %>%
      select(gene_id,GO_terms) %>%
      distinct() 
  ) %>%
  # filter(str_detect(GO_terms,"GO")) %>%
  # filter(!str_detect(gene_id,"Mac")) %>%
  mutate(GO_terms = str_replace_na(GO_terms,"NA")) %>%
  filter(GO_terms !='NA') %>%
  select(-gene_list) %>%
  distinct()


EV_MAB_protinfer_retrieved %>%
  filter(str_detect(gene_id,"Mac")) %>%
  tail()

# Global genes 

GO_terms_global_genes <-
  EV_mazia_specific %>%
  rbind(
    EV_bedadeti_specific,
    MA_specific,
    MB_specific
  ) %>%
  filter(str_detect(GO_terms,"GO")) %>%
    rbind(
      eggnog_MAB_retrieved,
      EV_MAB_protinfer_retrieved
    ) %>%
    distinct() 


# create an empty vector for saving output

# clustProfiler_GOenrich_out <- c()
# clustProfiler_GOenrich_out_allMapped <- c()
clustProfiler_GOenrich_output <- c()

# select target gene set
target_genes <-
  EV_bedadeti_specific %>%
  filter(str_detect(GO_terms,"GO"))  
 


# Extract global genes set
GO_terms_global_genes <-
  GO_terms_global_genes %>%
  filter(str_detect(GO_terms,"GO")) 

# List targe gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

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
  
  clustProfiler_GOenrich_output <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrich_output)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrich_output) = NULL

## plot 
# bedadeti_enrich_plot <-
clustProfiler_GOenrich_output %>%
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
  dplyr::filter(p.adjust < 0.05 ) %>%
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
  EV_mazia_specific %>%
  filter(str_detect(GO_terms,"GO")) #%>%

# Extract global genes set

GO_terms_global_genes <-
  GO_terms_global_genes %>%
  filter(str_detect(GO_terms,"GO")) 

# List target gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

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

  clustProfiler_GOenrich_output <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrich_output)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrich_output) = NULL

# plot
mazia_enrich_plot <-
clustProfiler_GOenrich_output %>%
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
  dplyr::filter(p.adjust < 0.05 ) %>%
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
  MA_specific %>%
  filter(str_detect(GO_terms,"GO")) #%>%

# Extract global genes set

GO_terms_global_genes <-
  GO_terms_global_genes %>%
  filter(str_detect(GO_terms,"GO")) 


# List target gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

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
  
  clustProfiler_GOenrich_output <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrich_output)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrich_output) = NULL


# plot 
# plot
MA_enrich_plot <-
clustProfiler_GOenrich_output %>%
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
  dplyr::filter(p.adjust < 0.05 ) %>%
  distinct()  %>%
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
  MB_specific %>%
  filter(str_detect(GO_terms,"GO")) #%>%

# Extract global genes set

GO_terms_global_genes <-
  GO_terms_global_genes %>%
  filter(str_detect(GO_terms,"GO")) 


# List target gene set
target_genes_list <-  unique(target_genes$gene_id)

enriched_out <-
  enricher(gene = target_genes_list, TERM2GENE = GO_terms_global_genes)

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
  
  clustProfiler_GOenrich_output <-
    enriched_out_df  %>% # add it to your list
    rbind(clustProfiler_GOenrich_output)    
  
}

# remove temporary files

rm(#GO_terms_global_genes,
  enriched_out,
  target_genes_list,
  enriched_out_df)

# tidyup output file
# row.names(clustProfiler_GOenrich_out) = NULL
row.names(clustProfiler_GOenrich_output) = NULL


# plot
MB_enrich_plot <-
clustProfiler_GOenrich_output %>%
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
  dplyr::filter(p.adjust < 0.05 ) %>%
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
clustProfiler_GOenrich_output %>%
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
clustProfiler_GOenrich_output %>%
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
  clustProfiler_GOenrich_output %>%
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
clustProfiler_GOenrich_output %>%
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

