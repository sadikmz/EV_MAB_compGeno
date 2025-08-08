#```{r}
#Set the working directory and load the input data
# setwd("/Volumes/Jamie_EXT/Research/Genomics/FolaProject/Maei/Analysis_V2_MaeiV4/Results/BlastSearches/6565-BigBlast_tblastn_16-12-2022/CandidateEffectorsFromHits")

## MA reaads 
panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,total)  %>% 
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,MAB_37,total)  %>%
  select(MAB_37,gene_id,total,genome) %>%
  filter(MAB_37 !=0) %>%

  # filter(MAB_37 ==1) %>%
  # filter(total < 65) %>%
  # filter(str_detect(gene_id,"prim"))
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>% 
  group_by(genome) %>%
  summarise(EV_proteins_with_no_MA_mapping_reads = n()) 
  
## MA reaads 
panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,total)  %>% 
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,total)  %>%
  # select(MAB_37,gene_id,total,genome) %>%
  filter(total >=1) %>%
  # filter(MAB_37 ==1) %>%
  # filter(total < 65) %>%
  # filter(str_detect(gene_id,"prim"))
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>% 
  group_by(genome) %>%
  summarise(EV_proteins_with_no_MB_mapping_reads = n()) 


## Common in all MA

panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,total)  %>% 
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,MAB_37,total)  %>%
  select(MAB_37,gene_id,total,genome) %>%
  filter(MAB_37 !=0) %>%
  
  # head()
  # mutate(
  #   total = MAB_69+MAB_68+MAB_67+MAB_66+MAB_65+MAB_64+MAB_63+MAB_62+MAB_61+MAB_60+MAB_59+MAB_58+
  #     MAB_57+MAB_56+MAB_55+MAB_54+MAB_53+MAB_52+MAB_51+MAB_50+MAB_49+MAB_48+MAB_47+MAB_46+MAB_45+
  #     MAB_44+MAB_43+MAB_42+MAB_41+MAB_40+MAB_39+MAB_38+MAB_36+MAB_35+MAB_34+MAB_33+MAB_32+MAB_31+
  #     MAB_30+MAB_29+MAB_28+MAB_27+MAB_26+MAB_25+MAB_24+MAB_23+MAB_22+MAB_21+MAB_20+MAB_19+MAB_18+
  #     MAB_17+MAB_16+MAB_15+MAB_14+MAB_13+MAB_12+MAB_11+MAB_10+MAB_9+MAB_8+MAB_7+MAB_6+MAB_5
  # ) %>%
  # filter(total ==64) %>%
  # head()
  # filter(MAB_37 ==1) %>%
  # filter(total < 65) %>%
  # filter(str_detect(gene_id,"prim"))
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>% 
  group_by(genome) %>%
  summarise(EV_proteins_with_no_MB_mapping_reads = n()) 


## common in all MB
panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,total)  %>% 
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-MAB_37,-total)  %>%
  # select(MAB_37,gene_id,total,genome) %>%
  mutate(
    total = MAB_69+MAB_68+MAB_67+MAB_66+MAB_65+MAB_64+MAB_63+MAB_62+MAB_61+MAB_60+MAB_59+MAB_58+
      MAB_57+MAB_56+MAB_55+MAB_54+MAB_53+MAB_52+MAB_51+MAB_50+MAB_49+MAB_48+MAB_47+MAB_46+MAB_45+
      MAB_44+MAB_43+MAB_42+MAB_41+MAB_40+MAB_39+MAB_38+MAB_36+MAB_35+MAB_34+MAB_33+MAB_32+MAB_31+
      MAB_30+MAB_29+MAB_28+MAB_27+MAB_26+MAB_25+MAB_24+MAB_23+MAB_22+MAB_21+MAB_20+MAB_19+MAB_18+
      MAB_17+MAB_16+MAB_15+MAB_14+MAB_13+MAB_12+MAB_11+MAB_10+MAB_9+MAB_8+MAB_7+MAB_6+MAB_5
  ) %>%
  filter(total ==64) %>%
  # head()
  # filter(MAB_37 ==1) %>%
  # filter(total < 65) %>%
  # filter(str_detect(gene_id,"prim"))
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>% 
  group_by(genome) %>%
  summarise(EV_proteins_with_no_MB_mapping_reads = n()) 


# AllData <-
panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                         TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  filter(total == 65) %>%
  # filter(str_detect(gene_id,"prim"))
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>%
  group_by(genome) %>%
  summarise(EV_proteins_with_no_MAB_mapping_reads = n()) %>%
  write.table(paste0(path_mapping,"EV_proteins_with_no_MAB_mapping_reads_fraction_coverage_less_0.25.txt"),
              col.names = T, row.names = F, quote = F)
  # sample_n(1000)


panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  filter(MAB_37 ==0) %>%
  filter(total >= 1) %>%
  # filter(str_detect(gene_id,"prim"))
  # select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,total)  %>%
  # select(MAB_37) %>%
  select(gene_id,genome) %>%
  # group_by(genome) %>%
  # summarise(EV_proteins_with_no_mapping_reads_only_from_MB = n()) %>%
  write.table(paste0(path_mapping,"EV_proteins_with_no_mapping_reads_only_from_MB_fraction_coverage_less_0.25.txt"),
              col.names = T, row.names = F, quote = F)
# # sample_n(1000)


AllData <-
panMAB_EV_MZBD_presence_absence_0_25.v2 %>%
  mutate(
    Predicted_Genes = case_when(str_detect(gene_id,"EVBD") ~ "E. ventricosum (Bedadeti)" ,
                                TRUE ~ "E. ventricosum (Mazia)")
  ) %>%
  # filter(total <= 30) %>%
  filter(total < 65) %>%
  # head(n=8) %>%
  
  # filter(str_detect(gene_id,"prim"))
  select(-contig_chr_name,-start,-end,-sum,-sum1,-sum2,-sum3,-total)
  

AllData %>%
  head()




# AllData %>%
#   nrow()
# read.csv(paste0(path_mapping,"AllCandidateEffectorSets_65Perc.clstr_heatmapdata.csv")) -> AllData 

#Select the data we require, remove the race and fsp. columns.
AllData %>%
  select(-genome,-Predicted_Genes) %>%
  column_to_rownames("gene_id") -> HeatmapData

# AllData %>%
#   select(-c("Global_Distribution","Lactucae_Distribution")) %>%
#   column_to_rownames("Cluster") -> HeatmapData

#Provide an overview of the data
dim(HeatmapData)

#Calculate the number of columns in the df.
HeatmapData_Columns <- ncol(HeatmapData)
#Convert to a binary data matrix
HeatmapData[,1:all_of(HeatmapData_Columns)]<-ifelse(HeatmapData[,1:all_of(HeatmapData_Columns)] >= 1, 1, 0) -> HeatmapData_Binary

#Order by col and row means
v_rowMeans <- rowMeans(HeatmapData_Binary)
v_colMeans <- colMeans(HeatmapData_Binary)
HeatmapData_Binary <- HeatmapData_Binary[order(v_rowMeans),order(v_colMeans)]

#```

#```{r}
#Clustering of input data
#Specifying method = binary as the dataset has been converted to a binary matrix.

#Clustering the Candidate effectors
HeatmapData_Binary %>%
  dist(, method = "binary") %>%
  hclust() -> CEC_Clustering
plot(CEC_Clustering)

#Clustering the fsp.
HeatmapData_Binary %>%
  t() %>% 
  dist(, method = "binary") %>%
  hclust() -> Fsp_Clustering
plot(Fsp_Clustering)

sessionInfo()

#```

#```{r}

#Seprate the distribution data from the matrix and create separate df. 

# AllData %>%
#   head()
  

AllData %>%
  select("gene_id", "Predicted_Genes") %>%
  column_to_rownames("gene_id") -> RowAnnotations #Add the distribution data

# AllData %>%
#   select(-genome) %>%
#   column_to_rownames("gene_id") -> HeatmapData

#Create a df for the isolates including race and fsp. 

AllData_colNames <- AllData %>%
  colnames() %>%
  str_replace("MAB_37","M. acuminata (2 samples)") %>%
  dput()

# AllData_colNames[1:65] 
data.frame(row.names = colnames(HeatmapData),
           Genomic_reads_source = AllData_colNames[1:65]) %>%
  mutate(Genomic_reads_source = case_when(str_detect(Genomic_reads_source,"MAB") ~ "M. balbisiana (63 samples)",
                                          TRUE ~ Genomic_reads_source)) -> ColAnnotations


#Create a list containing the colours we want to use for the annotations. 
ann_colours = list(
  Genomic_reads_source= c("M. balbisiana (63 samples)" = "orange2", "M. acuminata (2 samples)" = "blue"),
  # "Forma_specialis" = c(
  #   "apii" = "lightskyblue",
  #   "cepae" = "purple",
  #   "conglutinans" = "orange",
  #   "coriandrii" = "pink",
  #   "cubense" = "khaki1",
  #   "lactucae" = "green3",
  #   "lini" = "indianred",
  #   "lycopersici" = "red",
  #   "matthiolae" = "blue",
  #   "narcissus" = "yellow",
  #   "niveum" = "brown",
  #   "rapae" = "dark green",
  #   "rocket" = "light green",
  #   "vasinfectum" = "greenyellow",
  #   "non-pathogenic" = "black"),
  
  Predicted_Genes = c( "E. ventricosum (Bedadeti)" = "purple", "E. ventricosum (Mazia)" = "green"))

# Isolist = list("Fo47", "207A",	"Fus2",	"Fo5176",	"3-2",	"160527",	"UK0001",	"AJ516",	"AJ520",	"AJ592",	"AJ705",	"AJ856",	"AJ718",	"39",	"4287",	"AJ260",	"FON63",	"110407-3-1-1",	"Tf1208",	"TF1",	"AJ174")

RowAnnotations

ColAnnotations
#```

#```{r}

#Generate the heatmap. 
# PanMAB_read_EV_MZBD_prots_less_25 <- 
HeatmapData_Binary %>% 
  pheatmap(color = colorRampPalette(c("firebrick","grey90"))(2),
           legend = T,
           legend_breaks = c(0.75, 0.25),
           legend_labels = c("Absent","Present"),
           legend.cex = 25,
           drop_levels = F,
           cluster_rows = CEC_Clustering,
           cluster_cols = Fsp_Clustering,
           cellheight = 3,
           cellwidth = 50,
           annotation_row = RowAnnotations,
           annotation_col = ColAnnotations,
           annotation_colors = ann_colours,
           fontsize = 70,
           fontsize_row = fontsize,
           treeheight_row = 30,
           treeheight_col = 30,
           angle_col = 0,
           cutree_cols = 3,
           cutree_rows = 4,
           na_col = "white",
           border_color = "grey90",
           labels_col = NULL,
           show_colnames = FALSE,
           show_rownames = FALSE,
           filename = "./PanMAB_read_EV_MZBD_prots_less_25.4.jpeg",
           height = 50,
           width = 80)


# PanMAB_read_EV_MZBD_prots_less_25$gtable$grobs[[1]]$gp <- gpar(lwd=4)
# PanMAB_read_EV_MZBD_prots_less_25$gtable$grobs[[2]]$gp <- gpar(lwd=4)
# 
# png('PanMAB_read_EV_MZBD_prots_less_25.2.png', height = 50, width = 75)
# grid.newpage()
# grid.draw(PanMAB_read_EV_MZBD_prots_less_25$gtable)

HeatmapData_Binary %>% 
  pheatmap(color = colorRampPalette(c("firebrick","grey90"))(2),
           legend = T,
           legend_breaks = c(0.75, 0.25),
           legend_labels = c("Absent","Present"),
           drop_levels = F,
           cluster_rows = CEC_Clustering,
           cluster_cols = Fsp_Clustering,
           cellheight = 0.25,
           cellwidth = 10,
           annotation_row = RowAnnotations,
           annotation_col = ColAnnotations,
           annotation_colors = ann_colours,
           fontsize_row = 24,
           fontsize = 32,
           treeheight_row = 50,
           angle_col = 90,
           cutree_cols = 1,
           show_rownames = F,
           na_col = "white",
           border_color = "grey90",
           #labels_col = Isolist,
           filename = "PanMAB_read_EV_MZBD_prots.06.png",
           height = 500,
           width = 50)


citation("pheatmap")
#```

#Generate the heatmap. 
HeatmapData_Binary %>% 
  pheatmap(color = colorRampPalette(c("firebrick","grey90"))(2),
           legend = T,
           legend_breaks = c(0.75, 0.25),
           legend_labels = c("Absent","Present"),
           drop_levels = F,
           cluster_rows = CEC_Clustering,
           cluster_cols = Fsp_Clustering,
           cellheight = 6,
           cellwidth = 30,
           annotation_row = RowAnnotations,
           annotation_col = ColAnnotations,
           annotation_colors = ann_colours,
           fontsize = 25,
           fontsize_row = fontsize,
           treeheight_row = 30,
           treeheight_col = 30,
           angle_col = 90,
           cutree_cols = 2,
           show_rownames = F,
           na_col = "white",
           border_color = "grey90",
           # labels_col = FALSE,
           show_colnames = FALSE,
           filename = "./PanMAB_read_EV_MZBD_prots_less_25.png",
           height = 100,
           width = 100)


##
HeatmapData_Binary %>% 
  pheatmap(color = colorRampPalette(c("firebrick","grey90"))(2),
           legend = T,
           legend_breaks = c(0.75, 0.25),
           legend_labels = c("Absent","Present"),
           legend.cex = 4,
           drop_levels = F,
           cluster_rows = CEC_Clustering,
           cluster_cols = Fsp_Clustering,
           cellheight = 3,
           cellwidth = 40,
           annotation_row = RowAnnotations,
           annotation_col = ColAnnotations,
           annotation_colors = ann_colours,
           fontsize = 50,
           fontsize_row = fontsize,
           treeheight_row = 30,
           treeheight_col = 30,
           angle_col = 0,
           cutree_cols = 2,
           na_col = "white",
           border_color = "grey90",
           labels_col = NULL,
           show_colnames = FALSE,
           show_rownames = FALSE,
           filename = "./PanMAB_read_EV_MZBD_prots_less_25.3.jpeg",
           height = 50,
           width = 75)

