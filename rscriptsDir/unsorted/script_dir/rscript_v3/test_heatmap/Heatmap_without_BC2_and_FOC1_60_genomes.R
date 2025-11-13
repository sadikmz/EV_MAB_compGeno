---
  title: "Clustering of candidate effectors across all genomes"
output: html_notebook
---
  
  _Note:_ _if_ _you_ _wish_ _to_ _hide_ _the_ _code_, _please_ _select_ "Hide All Code" _the_ _dropdown_ _in_ _the_ _top_ _right_ _hand_ _corner._ 

---
  
  A set of candidate effectors identified using the pipeline developed in-house were blasted back (BLASTP, 1e-6 cutoff) against the 10 genome sequences used in the study, to indicate the presence and frequency of these candidate effectors across the genomes. 

---
  CURRENT WORKING DIRECTORY:
  ```{r}
#Jamie Pike 26/3/2021

#Pheatmap clusters by euclidean distance

##Used as a guide for this: https://www.youtube.com/watch?v=pTeTH9bz-_s  

#Set libraries 
#install.packages("tidyverse")
library(tidyverse)
library(pheatmap)
library(matrixStats)
library(knitr)
library(ggplot2)


getwd()

setwd("/Users/u1983390/Documents/Lab_Book/Candidate_Effectors/")

Effectors <- read.csv("transposed_list_of_effectors_in_genomes_Without_BC2_Foc_1_60.csv")
```
The table below indicates the number of BLASTP hits for each of the candidate effectors across all 10 genomes used in this study (Hits are counted where the sequences has => 905 identity). 
```{r}
# #Indicate the number of rows and columns in the csv. 
# print("nrow in csv:")
# nrow(Effectors) #Should equal the number of genomes.
# print("ncolv in csv:")
# ncol(Effectors) #Should equal the number of candidate effectors. 

#Set Heatmap data.
Effectors %>%
  select(1:11) -> selected_heatmap #n = 73 as we do not want to include the race data in the last column. 

print(selected_heatmap)
```

```{r} 
column_to_rownames(selected_heatmap, "Candidate_effectors") -> Heatmap_data
```
Initial Heatmap
===============
  
  The heatmap below is generated using the pheatmap package with the data from the BLASTP csv file indicated in the table above.
```{r fig1, fig.height = 20, fig.width = 7}
#Build Pheatmap

#(USE CMD + Shft + c to uncomment all text and look at plots which help to explain why the heatmap must be transformed.)
pheatmap(Heatmap_data) #note:
```
_Fig 1. clustering of candidate effectors across the Foc genomes._ 



There is no obvious clustering based on race, some TR4 isolates (e.g., FO_II5_V1) cluster with R1 isolates. This is due to the number of effectors found. 
This heatmap is hard to read due to the range of candidate effector numbers. E.g, Candid_Eff_67. Use plot below to visualize range in data.

Therefore, the data are log transformed, and now cluster more appropriately (race).

Final Heatmap
=============
  
  
  ```{r fig2, fig.height = 5, fig.width = 5 }
# #Therefore, the data are log transformed. 
# Effector_heatmap_data_log2_transformed <- log2(Heatmap_data+1) 
# Effector_heatmap_data_log2_transformed
# pheatmap(Effector_heatmap_data_log2_transformed) 
# # The genomes now cluster more appropriately. 
```
```{r fig3, fig.height = 3, fig.width = 8}
#Now to tidy up the heatmap.
#Generate even breaks in the scale bar for the heatmap.
#Symmetric_breaks <- seq(-max(abs(Effector_heatmap_data_log2_transformed)), max(abs(Effector_heatmap_data_log2_transformed)), length.out = 101)

Symmetric_breaks <- seq(-max(0), max(abs(Effector_heatmap_data_log2_transformed)), length.out = 101)

#Add race annotations 
data.frame(row.names = colnames(Heatmap_data),
           Race = c("TR4",	"R1",	"R4",	"R4",	"R1",	"TR4",	"Not identified",	"R1",	"TR4",	"TR4")) -> col_annotations
# print(col_annotations)

final_heatmap <- pheatmap(Effector_heatmap_data_log2_transformed,
                          #border_color = "black",
                          breaks = Symmetric_breaks,
                          annotation_col = col_annotations,
                          annotation_colors = list(
                            Race = c("Not identified" = "orange",
                                     "R1" = "white",
                                     "R4" = "grey60",
                                     "TR4" = "black")))

```
_Fig 2. clustering of candidate effectors across the Foc genomes. Data are log transformed and race indicated._

