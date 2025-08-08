AllData <- 
interproscan_GO_terms.v1.MAB 

AllData %>%
  head()

AllData %>%
  # select(-genome,-Predicted_Genes) %>%
  # mutate(description = str_to_sentence(description)) %>%
  column_to_rownames("description") -> HeatmapData


#```{r}

#Generate the heatmap. 

# heatmap_MAB <- 
AllData %>%
  filter(MA > 1) %>%
  mutate(description = str_to_sentence(description),
         description = str_replace(description,"dna","DNA"),
         description = str_replace(description,"Dna","DNA")) %>%
  column_to_rownames("description") %>%
  # head()
  pheatmap::pheatmap(#color = colorRampPalette(c("firebrick","grey90"))(2),
           # legend = T,
           # legend_breaks = c(0.75, 0.25),
           # legend_labels = c("Absent","Present"),
           legend.cex = 25,
           drop_levels = F,
           # cluster_rows = CEC_Clustering,
           # cluster_cols = Fsp_Clustering,
           # cellheight = 3,
           cellwidth = 50,
           # annotation_row = RowAnnotations,
           # annotation_col = ColAnnotations,
           # annotation_colors = ann_colours,
           fontsize = 13,
           # fontsize_row = fontsize,
           angle_col = 0,
           # cutree_cols = 3,
           # cutree_rows = 4,
           # na_col = "white",
           border_color = "grey90",
           # labels_col = NULL,
           show_colnames = TRUE,
           show_rownames = TRUE,
           # filename = "MAB_GO.heatmap.06.pdf",
           height = 6,
           width = 8,
           
           # hide dendrogram
           treeheight_row = 0,
           treeheight_col = 0
           )


# EV

AllData <- 
  interproscan_GO_terms.v1.EV

AllData %>%
  head()
# 
# AllData %>%
#   # select(-genome,-Predicted_Genes) %>%
#   # mutate(description = str_to_sentence(description)) %>%
#   column_to_rownames("description") -> HeatmapData


#```{r}

#Generate the heatmap. 

# heatmap_MAB <- 
AllData %>%
  dplyr::rename(Mazia = EV_mazia,
         Bedadeti = EV_bedadeti) %>%
  mutate(Mazia = str_replace_na(Mazia, "0"),
         Bedadeti = str_replace_na(Bedadeti, "0"),
         Mazia = as.numeric(Mazia),
         Bedadeti = as.numeric(Bedadeti)) %>% 
  filter(Mazia >= "0",
         Bedadeti > 1 |
           str_detect(description,"DNA") |
           description=="Biosynthetic process") %>%
  filter(!str_detect(description,"lysine Biosynthetic process via diaminopimelate"),
         !str_detect(description,"L-histidine Biosynthetic process"),
         !str_detect(description,"L-histidine Biosynthetic process")) %>%
  distinct() %>%
  select(description, Mazia, Bedadeti) %>%
  # filter(Mazia > 1) %>%
  # mutate(description = str_to_sentence(description),
  #        description = str_replace(description,"dna","DNA"),
  #        description = str_replace(description,"Dna","DNA")) %>%
  column_to_rownames("description") -> HeatmapData

HeatmapData <- HeatmapData[order(HeatmapData$Bedadeti),]

HeatmapData %>%
  # head()
  pheatmap::pheatmap(#color = colorRampPalette(c("firebrick","grey90"))(2),
    # legend = T,
    # legend_breaks = c(0.75, 0.25),
    # legend_labels = c("Absent","Present"),
    legend.cex = 25,
    drop_levels = F,
    # cluster_rows = CEC_Clustering,
    # cluster_cols = Fsp_Clustering,
    # cellheight = 3,
    cellwidth = 50,
    # annotation_row = RowAnnotations,
    # annotation_col = ColAnnotations,
    # annotation_colors = ann_colours,
    fontsize = 13,
    # fontsize_row = fontsize,
    angle_col = 0,
    # cutree_cols = 3,
    # cutree_rows = 4,
    # na_col = "white",
    border_color = "grey90",
    # labels_col = NULL,
    show_colnames = TRUE,
    show_rownames = TRUE,
    # filename = "EV_GO.heatmap.01.pdf",
    height = 6,
    width = 8,
    
    # hide dendrogram
    treeheight_row = 0,
    treeheight_col = 0
  )


## 
AllData %>%
  rename(Mazia = EV_mazia,
         Bedadeti = EV_bedadeti) %>%
  mutate(Mazia = str_replace_na(Mazia, "0"),
         Bedadeti = str_replace_na(Bedadeti, "0"),
         Mazia = as.numeric(Mazia),
         Bedadeti = as.numeric(Bedadeti)) %>% 
  filter(Bedadeti == 1,
           !str_detect(description,"DNA"),
           description !="Biosynthetic process",
         description !="respiratory electron transport chain") %>%
  filter(!str_detect(description,"lysine Biosynthetic process via diaminopimelate"),
         !str_detect(description,"L-histidine Biosynthetic process"),
         !str_detect(description,"L-histidine Biosynthetic process")) %>%
  distinct() %>%
  select(description, Mazia, Bedadeti) %>%
  # filter(Mazia > 1) %>%
  # mutate(description = str_to_sentence(description),
  #        description = str_replace(description,"dna","DNA"),
  #        description = str_replace(description,"Dna","DNA")) %>%
  column_to_rownames("description") -> HeatmapData

HeatmapData <- HeatmapData[order(HeatmapData$Bedadeti),]

HeatmapData %>%
  # head()
  pheatmap::pheatmap(#color = colorRampPalette(c("firebrick","grey90"))(2),
    # legend = T,
    # legend_breaks = c(0.75, 0.25),
    # legend_labels = c("Absent","Present"),
    legend.cex = 25,
    drop_levels = F,
    # cluster_rows = CEC_Clustering,
    # cluster_cols = Fsp_Clustering,
    # cellheight = 3,
    cellwidth = 50,
    # annotation_row = RowAnnotations,
    # annotation_col = ColAnnotations,
    # annotation_colors = ann_colours,
    fontsize = 13,
    # fontsize_row = fontsize,
    angle_col = 0,
    # cutree_cols = 3,
    # cutree_rows = 4,
    # na_col = "white",
    border_color = "grey90",
    # labels_col = NULL,
    show_colnames = TRUE,
    show_rownames = TRUE,
    # filename = "EV_GO.heatmap.00.pdf",
    height = 6,
    width = 8,
    
    # hide dendrogram
    treeheight_row = 0,
    treeheight_col = 0
  )

