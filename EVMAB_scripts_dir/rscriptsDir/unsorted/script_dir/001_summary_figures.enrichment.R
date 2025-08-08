# Plots for Overrepresentaiton analysis 


### Final plot
## Bedadeti
# bedadeti_enrich_plot <-
bedadeti_enrich_MABST <-
  clustProfiler_GOenrich_output_MAST %>%
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
  mutate(description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"carbohydrate metabolic process","Carbohydrate<br>metabolic process"),
         description = str_replace(description,"transmembrane transport","Transmembrane<br>transport"),
         description = str_replace(description,"proteolysis","Proteolysis"),
         description = str_replace(description,"translation","Translation")
  ) %>%
  dplyr::filter(aspects == 'Biological process') %>%
  dplyr::filter(p.adjust < 0.28 ) %>%
  distinct()#





bedadeti_enrich_MABST$description <- factor(bedadeti_enrich_MABST$description, 
                                            levels = c("Transmembrane<br>transport",
                                                       "regulation of transcription,<br>DNA-templated",
                                                       "Proteolysis", 
                                                       "Carbohydrate<br>metabolic process",
                                                       "signal transduction",
                                                       "Translation"
                                            ))


bedadeti_enrich_MABST_plot <-
  
  bedadeti_enrich %>%
  ggplot(aes(
    # x = reorder(description, p.adjust),
    x = description,
    
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
    axis.text.y = element_blank(),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "top",
    # legend.justification='left',
    legend.justification=c(0.93,0),
    axis.ticks.length.y = unit(.20, "cm")
    
    
  )

bedadeti_enrich_MABST_plot

# For editting plot

bedadeti_enrich_MABST_plot.v1 <-
  
  bedadeti_enrich %>%
  ggplot(aes(
    # x = reorder(description, p.adjust),
    x = description,
    
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
    axis.text.y = element_markdown(size=9, face="bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "top",
    # legend.justification='left',
    legend.justification=c(0.93,0),
    axis.ticks.length.y = unit(.20, "cm")
    
    
  )

bedadeti_enrich_MABST_plot.v1


# Mazia

mazia_enrich_MABST <-
  clustProfiler_GOenrich_output_MAST %>%
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
  # distinct(description)
  
  mutate(description = str_replace(description,"regulation of transcription, DNA-templated","regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"carbohydrate metabolic process","Carbohydrate<br>metabolic process"),
         description = str_replace(description,"transmembrane transport","Transmembrane<br>transport"),
         description = str_replace(description,"proteolysis","Proteolysis"),
         description = str_replace(description,"translation","Translation")
  ) %>%
  dplyr::filter(aspects == 'Biological process') %>%
  dplyr::filter(p.adjust < 0.2 ) #

mazia_enrich_MABST$description <- factor(mazia_enrich_MABST$description, 
                                         levels = c("Transmembrane<br>transport",
                                                    "regulation of transcription,<br>DNA-templated", 
                                                    "Proteolysis", 
                                                    "Carbohydrate<br>metabolic process", 
                                                    "DNA integration",
                                                    "Translation"))

mazia_enrich_MABST_plot <-
  mazia_enrich_MABST %>%
  ggplot(aes(
    # x = reorder(description, p.adjust),
    x = description,
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
  scale_x_discrete(position = "top") +
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
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "top",
    legend.justification=c(0.1,0),
    axis.ticks.length.y = unit(.20, "cm")
    
  )



mazia_enrich_MABST_plot

# for editting plot

mazia_enrich_MABST_plot.v1 <-
  mazia_enrich_MABST %>%
  ggplot(aes(
    # x = reorder(description, p.adjust),
    x = description,
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
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_blank(),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "top",
    legend.justification=c(0.1,0),
    axis.ticks.length.y = unit(.20, "cm")
    
  )


mazia_enrich_MABST_plot.v1


#MA specific 

MA_enrich_MABST <-
  clustProfiler_GOenrich_output_MAST %>%
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
  mutate(description = str_replace(description,"regulation of transcription, DNA-templated","Regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"recognition of pollen","Recognition of pollen"),
         description = str_replace(description,"protein phosphorylation","Protein phosphorylation"),
         description = str_replace(description,"intracellular protein transport","Intracellular<br>protein transport"),
         description = str_replace(description,"defense response","Defense response")
  ) %>%
  # dplyr::filter(aspects == 'Biological process') %>%
  dplyr::filter(p.adjust < 0.05 ) %>% 
  dplyr::filter(aspects == 'Biological process') %>%
  distinct()



MA_enrich_MABST$description <- factor(MA_enrich_MABST$description, 
                                      levels = c("Protein phosphorylation", 
                                                 "Regulation of transcription,<br>DNA-templated", 
                                                 "Defense response", 
                                                 "Intracellular<br>protein transport", 
                                                 "Recognition of pollen"))

MA_enrich_MABST_plot <-
  MA_enrich_MABST %>%
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
  scale_x_discrete(position = "top") +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    # axis.text.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    legend.justification=c(0.1,0),
    axis.ticks.length.y = unit(.20, "cm")
  )

MA_enrich_MABST_plot

# for editting plot 

MA_enrich_MABST_plot.v1 <- 
  MA_enrich_MABST %>%
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
  scale_x_discrete(position = "top") +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    # axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    legend.justification=c(0.1,0),
    axis.ticks.length.y = unit(.20, "cm")
  )


#MB specific 

MB_enrich_MBAST <-
  clustProfiler_GOenrich_output_MAST %>%
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
  dplyr::filter(description != "NA") %>%
  mutate(description = str_replace(description,"regulation of transcription, DNA-templated","Regulation of transcription,<br>DNA-templated"),
         description = str_replace(description,"transmembrane transport","Transmembrane transport"),
         description = str_replace(description,"protein phosphorylation","Protein phosphorylation"),
         description = str_replace(description,"translation","Translation"),
         description = str_replace(description,"defense response","Defense response")
  ) %>%
  dplyr::filter(p.adjust < 0.05 ) %>% 
  dplyr::filter(aspects == 'Biological process') %>%
  distinct()



MB_enrich_MBAST$description <- factor(MB_enrich_MBAST$description, 
                                      levels = c("Protein phosphorylation", 
                                                 "Regulation of transcription,<br>DNA-templated", 
                                                 "Defense response", 
                                                 "Translation", 
                                                 "Transmembrane transport"))

MB_enrich_MBAST_plot <-
  MB_enrich_MBAST %>%
  ggplot(aes(
    x = description,
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
  # scale_x_discrete(position = "top") +
  
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
    axis.text.y = element_blank(),
    # axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    # legend.justification='right'
    legend.justification=c(.9,0),
    axis.ticks.length.y = unit(.20, "cm")
  )




MB_enrich_MBAST_plot

# For editting plot 

MB_enrich_MBAST_plot.v1 <-
  MB_enrich_MBAST %>%
  ggplot(aes(
    x = description,
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
  # scale_x_discrete(position = "top") +
  
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
    # axis.text.y = element_blank(),
    axis.text.y = element_markdown(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    # legend.justification='right'
    legend.justification=c(.9,0),
    axis.ticks.length.y = unit(.20, "cm")
  )

library(patchwork)


#MABST
mazia_enrich_MABST_plot | bedadeti_enrich_MABST_plot
ggsave(paste0(DATA_DIR, "plot_out/EV_specific.enrichallMapped.MABST.synteny.pdf"), 
       width=8, height=5)

mazia_enrich_MABST_plot.v1 | bedadeti_enrich_MABST_plot.v1
ggsave(paste0(DATA_DIR, "plot_out/EV_specific.enrichallMapped.MABST.synteny.01.pdf"), 
       width=8, height=5)


MA_enrich_MABST_plot|MB_enrich_MBAST_plot
ggsave(paste0(DATA_DIR, "plot_out/MAB_specific.enrichallMapped.MABST.synteny.pdf"), 
       width=8, height=5)

MA_enrich_MABST_plot.v1|MB_enrich_MBAST_plot.v1
ggsave(paste0(DATA_DIR, "plot_out/MAB_specific.enrichallMapped.MABST.synteny.1.pdf"), 
       width=8, height=5)

ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.02.pdf"), width=11, height=4)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.02.png"), width=11, height=4)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.03.pdf"), width=11, height=4.5)
ggsave(paste0(DATA_DIR, "plot_out/EV_MAB.enrichallMapped.synteny.03.png"), width=11, height=4.5)