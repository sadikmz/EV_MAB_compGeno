# BiocManager::install("clusterProfiler")
library(tidyverse)
library(clusterProfiler)
library(ggtext)
library(rcartocolor)


# Set of identified genes specific to each species/genotype 
EV_MAB_custered_genes.v5 %>%
  head()

# Join with GO terms annotation from interproscan

EV_MAB_custered_genes.v5 %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      distinct(gene_id, GO_terms)
  ) %>%
  head()

interproscan_GO_terms <-
  EV_MAB_custered_genes.v5 %>%
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


bedadeti_specific_plot <-
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "EV_bedadeti_specific") %>%
  mutate(
    description = str_replace(description, "isopentenyl diphosphate biosynthetic process...", "isopentenyl diphosphate biosynthetic process"),
    # description = str_replace(description, "transmembrane transport", "transmembrane stransport")
  ) %>%
  # head()
  filter(!str_detect(description, "proteolysis"),
           description != "lysine biosynthetic process via diaminopimelate", 
         description != "iron-sulfur cluster assembly",
         description != "iron-sulfur cluster assembly",
         description != "transsulfuration",
         description != "amino acid metabolic process",
         description != "L-histidine biosynthetic process") %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(Count >=  2 ) %>%
  # distinct(description)

  rbind(
    interproscan_GO_terms.v1 %>% 
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "EV_bedadeti_specific") %>%
      dplyr::filter(description != "NA") %>%
      dplyr::filter(Count <=  1 ) %>%
      filter(description == "DNA transposition" )
  ) %>%
  mutate(
    description = str_replace(description, "transmembrane transport","Transmembrane transport"),
    description = str_replace(description, "regulation of DNA-templated transcription","Regulation of DNA-templated transcription"),
    description = str_replace(description, "phosphorelay signal transduction system","Phosphorelay signal transduction system"),
    description = str_replace(description, "phosphorylation","Phosphorylation"),
    description = str_replace(description, "carbohydrate metabolic process","Carbohydrate metabolic process"),
    description = str_replace(description, "amino acid metabolic process","Amino acid metabolic process"),
    description = str_replace(description, "glycolytic process","Glycolytic process"),
    description = str_replace(description, "gluconeogenesis","Gluconeogenesis"),
    description = str_replace(description, "glutathione metabolic process","Glutathione metabolic process"),
    description = str_replace(description, "transsulfuration","Transsulfuration"),
    description = str_replace(description, "protein folding","Protein folding"),
    
  ) %>%

  ggplot(aes(x = reorder(description, Count),y = Count)) +
  geom_col(width = 0.8,
           fill = "blue",
           alpha = 0.25
           ) +
  scale_y_continuous(breaks =c(1,3,5,7,8))+
  theme_classic() +
  # theme_classic()+
  coord_flip() +
  labs(
    y = "Count of genes" ,
    x = "GO terms description for gene set",
  ) +
  # scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )


interproscan_GO_terms %>%
  head()


# Mazia
mazia_specific_plot <- 
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "EV_mazia_specific") %>%
  filter(!str_detect(description, "respiratory electron transport chain"),
         !str_detect(description, "phosphoenolpyruvate"),
         !str_detect(description, "photosynthesis"),
         !str_detect(description, "ATP metabolic process"),
         !str_detect(description, "cellular oxidant detoxification"),
         !str_detect(description, "microtubule-based movement"),
         !str_detect(description, "photosynthetic electron transport chain"),
         !str_detect(description, "vesicle-mediated transport"),
         !str_detect(description, "actin filament organization"),
         !str_detect(description, "nucleosome assembly"),
         !str_detect(description, "regulation of proton transport"),
         !str_detect(description, "translation")
         ) %>%
  mutate(
    description = str_replace(description, "intracellular signal transduction","Intracellular signal transduction"),
    description = str_replace(description, "regulation of DNA-templated transcription","Regulation of DNA-templated transcription"),
    description = str_replace(description, "signal transduction system","Signal transduction system"),
    description = str_replace(description, "carbohydrate metabolic process","Carbohydrate metabolic process"),
    description = str_replace(description, "proteasome assembly","Proteasome assembly"),
    description = str_replace(description, "glycolytic process","Glycolytic process"),
    description = str_replace(description, "proton transmembrane transport","Proton transmembrane transport"),
    description = str_replace(description, "lipid metabolic process","Lipid metabolic process"),
    description = str_replace(description, "biosynthetic process","Biosynthetic process"),
    description = str_replace(description, "protein phosphorylation","Protein phosphorylation")
  ) %>%
  ggplot(aes(x = reorder(description, Count),y = Count)) +
  geom_col(width = 0.8,
           fill = "blue",
           alpha = 0.25
  ) +
  scale_y_continuous(breaks =c(1,2))+
  theme_classic() +
  coord_flip() +
  labs(
    y = "Count of genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )


#MA

MA_specific_plot <- 
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "MA_specific") %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(Count >  1 ) %>%
  filter(!str_detect(description, "respiratory electron transport chain"),
         !str_detect(description, "phosphoenolpyruvate"),
         !str_detect(description, "photosynthesis"),
         !str_detect(description, "ATP metabolic process"),
         !str_detect(description, "cellular oxidant detoxification"),
         !str_detect(description, "microtubule-based movement"),
         !str_detect(description, "photosynthetic electron transport chain"),
         !str_detect(description, "vesicle-mediated transport"),
         !str_detect(description, "actin filament organization"),
         !str_detect(description, "nucleosome assembly"),
         !str_detect(description, "regulation of proton transport"),

         !str_detect(description, "translation"),
         # MA addtional
         !str_detect(description, "vacuolar transport"),
         !str_detect(description, "folding"),
         !str_detect(description, "prolyl-tRNA aminoacylation"),
         
         !str_detect(description, "RNA processing"),
         !str_detect(description, "protein ubiquitination"),
         !str_detect(description, "negative regulation of organ growth"),
         !str_detect(description, "mitochondrial")
  ) %>%
  rbind(
    interproscan_GO_terms.v1 %>% 
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "MA_specific") %>%
      dplyr::filter(description != "NA") %>%
  dplyr::filter(Count == 1 ) %>%
  filter(str_detect(description, "carbohydrate metabolic process")|
         str_detect(description, "vernalization response"))) %>%
  ggplot(aes(x = reorder(description, Count),y = Count)) +
  geom_col(width = 0.8,
           fill = "red",
           alpha = 0.25) +
  scale_y_continuous(breaks =c(1,5,10,15,20,25))+
  theme_classic() +
  coord_flip() +
  labs(
    y = "Count of genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )


# MB

MB_specific_plot <-
  interproscan_GO_terms.v1 %>%
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "MB_specific") %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(Count >  2 ) %>%
  filter(!str_detect(description, "respiratory electron transport chain"),
         !str_detect(description, "phosphoenolpyruvate"),
         !str_detect(description, "photosynthesis"),
         !str_detect(description, "ATP metabolic process"),
         !str_detect(description, "cellular oxidant detoxification"),
         !str_detect(description, "microtubule-based movement"),
         !str_detect(description, "photosynthetic electron transport chain"),
         !str_detect(description, "vesicle-mediated transport"),
         !str_detect(description, "actin filament organization"),
         !str_detect(description, "nucleosome assembly"),
         !str_detect(description, "regulation of proton transport"),

         !str_detect(description, "translation"),
         # MA addtional
         !str_detect(description, "vacuolar transport"),
         !str_detect(description, "RNA processing"),
         !str_detect(description, "protein ubiquitination"),
         !str_detect(description, "negative regulation of organ growth"),
         !str_detect(description, "mitochondrial")) %>%
  rbind(
    interproscan_GO_terms.v1 %>%
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "MB_specific") %>%
      dplyr::filter(description != "NA") %>%
      dplyr::filter(Count <=  2 ) %>%

      filter(str_detect(description, "carbohydrate metabolic process")|
               str_detect(description, "DNA-templated transcription")|
               str_detect(description, "tricarboxylic acid cycle")|
               str_detect(description, "vernalization response"))) %>%
    
  ggplot(aes(x = reorder(description, Count),y = Count)) +
  # geom_col(width = 0.8, color = "white", alpha = 0.1) +
  geom_col(width = 0.8,
           fill = "red",
           alpha = 0.25) +
  theme_classic() +
  coord_flip() +
  labs(
    y = "Count of genes" ,
    x = "GO terms description for gene set (Specific to MABS specific)",
  ) +
  # scale_fill_gradientn(colours = terrain.colors(10), name  = "log10 pvalue") +
  scale_y_continuous(breaks =c(1,3,5,7,9,12))+
      
  theme(
    axis.title.x = element_markdown(size = 9, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_markdown(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )

# combine plots
library(patchwork)
mazia_specific_plot+bedadeti_specific_plot
MA_specific_plot+MB_specific_plot
