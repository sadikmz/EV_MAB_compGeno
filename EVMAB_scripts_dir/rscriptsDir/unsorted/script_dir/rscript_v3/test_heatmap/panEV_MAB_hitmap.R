#```{r}
#Set the working directory and load the input data
# setwd("/Volumes/Jamie_EXT/Research/Genomics/FolaProject/Maei/Analysis_V2_MaeiV4/Results/BlastSearches/6565-BigBlast_tblastn_16-12-2022/CandidateEffectorsFromHits")

read.csv(paste0(DATA_DIR,"../script_dir/rscript_v3/test_heatmap/AllCandidateEffectorSets_65Perc.clstr_heatmapdata.csv")) -> AllData 

read.csv(paste0(DATA_DIR,"../script_dir/rscript_v3/test_heatmap/AllCandidateEffectorSets_65Perc.clstr_heatmapdata.csv")) -> AllData 

list.files(paste0(DATA_DIR,"../script_dir/rscript_v3/test_heatmap/"))

#Select the data we require, remove the race and fsp. columns.
AllData %>% 
  select(-c("Global_Distribution","Lactucae_Distribution")) %>%
  column_to_rownames("Cluster") -> HeatmapData

AllData %>%
  select(-c("Global_Distribution","Lactucae_Distribution")) %>%
  column_to_rownames("Cluster") %>%
  head()

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
AllData %>% 
  select("Cluster", "Lactucae_Distribution") %>%
  column_to_rownames("Cluster") -> RowAnnotations #Add the distribution data


RowAnnotations %>%
  head()

#Create a df for the isolates including race and fsp. 
data.frame(row.names = colnames(HeatmapData),
           Race = c(NA, NA, NA, NA, NA, NA, NA,  "Race 4", "Race 1", "Race 4", "Race 4", "Race 1", "Race 1", NA, NA, NA, NA, NA, NA, NA,NA),
           "Forma_specialis"  = c("non-pathogenic", "apii", "cepae", "conglutinans", "coriandrii", "cubense", "cubense",  "lactucae", "lactucae", "lactucae", "lactucae", "lactucae", "lactucae", "lini",
                                  "lycopersici", "matthiolae", "narcissus", "niveum", "rapae", "vasinfectum",
                                  "rocket")) -> ColAnnotations

ColAnnotations %>%
  head()

#Create a list containing the colours we want to use for the annotations. 
ann_colours = list(
  Race= c("Race 1" = "orange2", "Race 4" = "blue"),
  "Forma_specialis" = c(
    "apii" = "lightskyblue",
    "cepae" = "purple",
    "conglutinans" = "orange",
    "coriandrii" = "pink",
    "cubense" = "khaki1",
    "lactucae" = "green3",
    "lini" = "indianred",
    "lycopersici" = "red",
    "matthiolae" = "blue",
    "narcissus" = "yellow",
    "niveum" = "brown",
    "rapae" = "dark green",
    "rocket" = "light green",
    "vasinfectum" = "greenyellow",
    "non-pathogenic" = "black"), 
  
  Lactucae_Distribution = c( "Race 1 Specific" = "tan1", "Race 4 Specific" = "skyblue3", "Shared" = "green", "Not Race Specific" = "grey"))

Isolist = list("Fo47", "207A",	"Fus2",	"Fo5176",	"3-2",	"160527",	"UK0001",	"AJ516",	"AJ520",	"AJ592",	"AJ705",	"AJ856",	"AJ718",	"39",	"4287",	"AJ260",	"FON63",	"110407-3-1-1",	"Tf1208",	"TF1",	"AJ174")

RowAnnotations

ColAnnotations
#```

#```{r}

#Generate the heatmap. 
HeatmapData_Binary %>%
  pheatmap(color = colorRampPalette(c("grey90", "firebrick"))(2),
           legend = T,
           legend_breaks = c(0.75, 0.25),
           legend_labels = c("Present", "Absent"),
           drop_levels = F,
           cluster_rows = CEC_Clustering,
           cluster_cols = Fsp_Clustering,
           cellheight = 5.5,
           cellwidth = 60,
           annotation_row = RowAnnotations,
           annotation_col = ColAnnotations,
           annotation_colors = ann_colours,
           fontsize_row = 2.5,
           fontsize = 24,
           treeheight_row = 10,
           angle_col = 90,
           cutree_cols = 14,
           show_rownames = F,
           na_col = "white",
           border_color = "grey90",
           #labels_col = Isolist,
           filename = "./FolaEffectorsHeatmap-Draft2.pdf",
           height = 30,
           width = 25)

citation("pheatmap")
#```


