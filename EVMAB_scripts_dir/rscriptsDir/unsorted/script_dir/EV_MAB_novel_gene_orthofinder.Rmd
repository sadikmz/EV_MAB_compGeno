```{r libraries}
library(tidyverse)
library(cowplot)
library(ggVennDiagram)
library(tidyverse)
library(glue)
library(ggtext)
library(ggrepel)
```

# Read input data
## Assigned and unassigne Orthogroups
```{r}
## assigend and unassinged orthogroups or EV, EG and Musa genomes 


EV_AED20_MAB.OG_GeneCount <-
  read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups//Orthogroups.GeneCount.tsv", sep = '\t', header = T) %>%  
  mutate(
    Ensete_ventricosum_mazia = case_when(Ensete_ventricosum_mazia >=1 ~ 1,
                                         TRUE ~ 0),
    Ensete_ventricosum_bedadeti = case_when(Ensete_ventricosum_bedadeti >=1 ~ 1,
                                            TRUE ~ 0),
    Musa_acuminata = case_when(Musa_acuminata >=1 ~ 1,
                               TRUE ~ 0),
    Musa_balbisiana = case_when(Musa_balbisiana >=1 ~ 1,
                                TRUE ~ 0)
    # Musa_schizocarpa = case_when(Musa_schizocarpa >=1 ~ 1,
    #                              TRUE ~ 0)
  ) %>%
  
  ## add unassigned OGs 
  bind_rows(
    read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/Orthogroups_UnassignedGenes.tsv", header = T, sep = '\t') %>%
      mutate(
        Ensete_ventricosum_mazia = case_when(str_detect(Ensete_ventricosum_mazia,"\\w+") ~ 1,
                                             TRUE ~ 0),
        Ensete_ventricosum_bedadeti = case_when(str_detect(Ensete_ventricosum_bedadeti, "\\w+") ~ 1,
                                                TRUE ~ 0),
        Musa_acuminata = case_when(str_detect(Musa_acuminata, "\\w+") ~ 1,
                                   TRUE ~ 0),
        Musa_balbisiana = case_when(str_detect(Musa_balbisiana, "\\w+") ~ 1,
                                    TRUE ~ 0)
        # Musa_schizocarpa = case_when(str_detect(Musa_schizocarpa, "\\w+") ~ 1,
        #                              TRUE ~ 0)
      )) %>% select(-Total)


```  



## Prepare dat for Ven Diagram
```{r}
dir.create("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/")


EV_AED20_MAB.OG_GeneCount %>%
  filter(Ensete_ventricosum_mazia == 1) %>% 
  select(Orthogroup) %>%
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Ensete_ventricosum_mazia.OG.txt", col.names = F, sep = '\t', quote = F, row.names = F)

EV_AED20_MAB.OG_GeneCount %>%
  filter(Ensete_ventricosum_bedadeti == 1) %>% 
  select(Orthogroup) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Ensete_ventricosum_bedadeti.OG.txt", col.names = F, sep = '\t', quote = F, row.names = F)


EV_AED20_MAB.OG_GeneCount %>%
  filter(Musa_acuminata == 1) %>% 
  select(Orthogroup) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Musa_acuminata.OG.txt", col.names = F, sep = '\t', quote = F, row.names = F)

EV_AED20_MAB.OG_GeneCount %>%
  filter(Musa_balbisiana == 1) %>% 
  select(Orthogroup) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Musa_balbisiana.OG.txt", col.names = F, sep = '\t', quote = F, row.names = F)



# Figure of main venn diagrams

OG_Species = c("Ensete_ventricosum_bedadeti","Ensete_ventricosum_mazia","Musa_acuminata","Musa_balbisiana")

OG_list_VennDiag_EV_EG_MUSA <- list()

OG_list_VennDiag_EV_EG_MUSA[["Ensete ventricosum\n(Bedadeti)"]] <- read.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Ensete_ventricosum_bedadeti.OG.txt", header = F, stringsAsFactors = FALSE)$V1
OG_list_VennDiag_EV_EG_MUSA[["Ensete ventricosum\n(Mazia)"]] <- read.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Ensete_ventricosum_mazia.OG.txt", header = F, stringsAsFactors = FALSE)$V1

OG_list_VennDiag_EV_EG_MUSA[["Musa acuminata"]] <- read.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Musa_acuminata.OG.txt", header = F, stringsAsFactors = FALSE)$V1
OG_list_VennDiag_EV_EG_MUSA[["Musa balbisiana"]] <- read.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Musa_balbisiana.OG.txt", header = F, stringsAsFactors = FALSE)$V1

```

# Plot Ven diagra 
```{r vendiag}
library(VennDiagram)

ggVennDiagram(x = OG_list_VennDiag_EV_EG_MUSA, 
              size = 4.3) +
  scale_fill_gradient(low="light grey", high = "red", limits = c(0, 18000)) +
  labs(fill="Count") +
  # ggtitle("Shared Orthogroups (OGs)") +
  theme(plot.title = element_text(hjust = 0.6),
        legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/EV_AED20_MAB_VenDiag.tiff", width=8, height=4)
```



## summary statistics 
```{r summary_stat}
### statistics 
read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Comparative_Genomics_Statistics/Statistics_Overall.tsv", sep = '\t', header = F) %>%
  slice(1:18) %>% 
  mutate(V1 = str_replace_all(V1, "Number", "#"),
         V1 = str_replace(V1, "Percentage", "%"),
         V1 = str_replace(V1, "orthogroup", "OG")) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Comparative_Genomics_Statistics/Statistics_Overall.txt", sep = '\t', col.names = F, row.names = F, quote = F)
```

# Prepare data for duplication in species species OGs 
```{r dups_species_spec}
# Genome specific OGs
species_specific_OGs_all <-
  EV_AED20_MAB.OG_GeneCount %>% 
  filter(
    Ensete_ventricosum_bedadeti == 1,
    Ensete_ventricosum_mazia==0, 
    Musa_acuminata == 0,
    Musa_balbisiana ==0
  ) %>%
  select(Orthogroup) %>% 
  mutate(genome="EV_bedadeti") %>% 
  bind_rows(
    EV_AED20_MAB.OG_GeneCount %>% 
      filter(
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia==1, 
        Musa_acuminata == 0,
        Musa_balbisiana ==0
        
      ) %>%
      select(Orthogroup) %>% 
      mutate(genome="EV_mazia"),
    EV_AED20_MAB.OG_GeneCount %>% 
      filter(
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia==0, 
        Musa_acuminata == 1,
        Musa_balbisiana ==0
        
      ) %>%
      select(Orthogroup) %>% 
      mutate(genome="MA"),
    EV_AED20_MAB.OG_GeneCount %>% 
      filter(
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia==0, 
        Musa_acuminata == 0,
        Musa_balbisiana ==1
        
      ) %>%
      select(Orthogroup) %>% 
      mutate(genome="MB"))


duplications_all <-
  read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Gene_Duplication_Events/Duplications.tsv") 
```

### plot
```{r plot_dup_species_specific}
# duplications_all %>% 
species_specific_OGs_all %>% 
  # filter(genome == "EV_bedadeti") %>% 
  left_join(
    duplications_all
  ) %>% distinct() %>% 
  mutate(Species.Tree.Node = str_replace_na(Species.Tree.Node,"NA")) %>%
  filter(Species.Tree.Node != "NA") %>%
  select(Orthogroup, genome) %>% 
  distinct() %>% 
  group_by(genome) %>%
  summarise( duplications_OG = n()) %>%
  left_join(
    species_specific_OGs_all %>%
      group_by(genome) %>%
      summarise(total = n()) %>% 
      arrange(total) 
  ) %>%
  pivot_longer(!genome, names_to = "OG_category", values_to = "value") %>% 
  mutate(OG_category = str_replace(OG_category, "duplications_OG", "OGs with duplications event"),
         OG_category = str_replace(OG_category, "total", "Genome specific OGs")) %>% 
  # rename('OGs category' = OG_category) %>%
  ggplot(aes(x= factor(genome, levels = c("EV_bedadeti","EV_mazia", "MB","MA")),
             y=value, fill= OG_category ))+
  geom_col( position = "dodge" )+
  scale_y_continuous(breaks =c(1139, 1927, 2477,6095))+
  scale_x_discrete(breaks =c( "MB", "MA", "EV_bedadeti", "EV_mazia"),
                   labels = c("MB","MA",  "EV (Bedadeti)","EV (Mazia)"))+
  
  coord_flip()+
  # theme_bw() +
  labs(
    y =" Number of genome specifc orthogroups (OGs) or<br> OGs that went one or more duplications event",
    # x = "Genomes"
  )+
  theme(
    axis.title.x = element_markdown(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_markdown(size = 12),
    # axis.title.x
    panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 12),
    legend.position = c(0.7,0.9)
  )
```


# Filter genes in duplications event 
```{r dup_speceis_specific_genes}
## bedadeti

species_specific_OGs_all %>% 
  # filter(genome == "EV_bedadeti") %>% 
  left_join(
    duplications_all
  ) %>% 
  mutate(Species.Tree.Node = str_replace_na(Species.Tree.Node,"NA")) %>%
  filter(Species.Tree.Node != "NA") %>% 
  distinct() %>% 
  filter(genome == "EV_bedadeti" ) %>% 
  select(Orthogroup, Genes.1)  %>% 
  mutate(Genes.1 = str_remove_all(Genes.1,"Ensete_ventricosum_bedadeti_"),
         Genes.1 = str_remove_all(Genes.1, "\\s\\w+\\s+\\w+_0.\\d+\\s+\\w+_0.\\d+"),
         Genes.1 = str_replace_all(Genes.1,"\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""),
         # ,
         # Genes.1 = str_replace_all(Genes.1,"\\s+\\d+|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""
         
         
         Genes.1 = str_replace_all(Genes.1,"\\d+\\|\\d+\\|\\d+\\|", ""),
         Genes.1 = str_replace_all(Genes.1,"_\\d+", "_"),
         Genes.1 = str_replace_all(Genes.1,"_.\\d+", "_"),
         Genes.1 = str_replace_all(Genes.1,"QI_", ""),
         Genes.1 = str_replace_all(Genes.1,"\\|\\d+", ""),
         Genes.1 = str_replace_all(Genes.1,"\\|--\\d+", "")
  ) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Gene_Duplication_Events/bedadeti_unique_dup_gene1.txt", col.names = F, row.names = F, quote = F, sep='\t')
# distinct()  %>% 
# mutate(genome="EV_beadeti") %>% 

species_specific_OGs_all %>% 
  filter(genome == "EV_bedadeti") %>% 
  left_join(
    duplications_all
  ) %>% distinct() %>% 
  mutate(Species.Tree.Node = str_replace_na(Species.Tree.Node,"NA")) %>%
  filter(Species.Tree.Node != "NA") %>% 
  select(Orthogroup, Genes.2)  %>% 
  mutate(Genes.2 = str_remove_all(Genes.2,"Ensete_ventricosum_bedadeti_"),
         Genes.2 = str_remove_all(Genes.2, "\\s\\w+\\s+\\w+_0.\\d+\\s+\\w+_0.\\d+"),
         Genes.2 = str_replace_all(Genes.2,"\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""),
         # ,
         # Genes.1 = str_replace_all(Genes.1,"\\s+\\d+|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""
         Genes.2 = str_replace_all(Genes.2,"\\d+\\|\\d+\\|\\d+\\|", ""),
         Genes.2 = str_replace_all(Genes.2,"_\\d+", "_"),
         Genes.2 = str_replace_all(Genes.2,"_.\\d+", "_"),
         Genes.2 = str_replace_all(Genes.2,"QI_", ""),
         Genes.2 = str_replace_all(Genes.2,"\\|\\d+", ""),
         Genes.2 = str_replace_all(Genes.2,"\\|--\\d+", "")
  ) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Gene_Duplication_Events/bedadeti_unique_dup_gene2.txt", col.names = F, row.names = F, quote = F, sep = '\t')
# distinct()  %>% 
# mutate(genome="EV_beadeti") %>% 

## Mazia 
species_specific_OGs_all %>% 
  # filter(genome == "EV_bedadeti") %>% 
  left_join(
    duplications_all
  ) %>% 
  mutate(Species.Tree.Node = str_replace_na(Species.Tree.Node,"NA")) %>%
  filter(Species.Tree.Node != "NA") %>% 
  distinct() %>% 
  filter(genome == "EV_mazia" ) %>% 
  select(Orthogroup, Genes.1)  %>% 
  mutate(Genes.1 = str_remove_all(Genes.1,"Ensete_ventricosum_mazia_"),
         Genes.1 = str_remove_all(Genes.1, "\\s\\w+\\s+\\w+_0.\\d+\\s+\\w+_0.\\d+"),
         Genes.1 = str_replace_all(Genes.1,"\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""),
         # ,
         # Genes.1 = str_replace_all(Genes.1,"\\s+\\d+|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""
         
         
         Genes.1 = str_replace_all(Genes.1,"\\d+\\|\\d+\\|\\d+\\|", ""),
         Genes.1 = str_replace_all(Genes.1,"_\\d+", "_"),
         Genes.1 = str_replace_all(Genes.1,"_.\\d+", "_"),
         Genes.1 = str_replace_all(Genes.1,"QI_", ""),
         Genes.1 = str_replace_all(Genes.1,"\\|\\d+", ""),
         Genes.1 = str_replace_all(Genes.1,"\\|--\\d+", "")
  ) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Gene_Duplication_Events/mazia_unique_dup_gene1.txt", col.names = F, row.names = F, quote = F, sep = '\t')
# distinct()  %>% 
# mutate(genome="EV_beadeti") %>% 

species_specific_OGs_all %>% 
  filter(genome == "EV_mazia") %>% 
  left_join(
    duplications_all
  ) %>% distinct() %>% 
  mutate(Species.Tree.Node = str_replace_na(Species.Tree.Node,"NA")) %>%
  filter(Species.Tree.Node != "NA") %>% 
  select(Orthogroup, Genes.2)  %>% 
  mutate(Genes.2 = str_remove_all(Genes.2,"Ensete_ventricosum_mazia_"),
         Genes.2 = str_remove_all(Genes.2, "\\s\\w+\\s+\\w+_0.\\d+\\s+\\w+_0.\\d+"),
         Genes.2 = str_replace_all(Genes.2,"\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""),
         # ,
         # Genes.1 = str_replace_all(Genes.1,"\\s+\\d+|\\d+\\|\\d+\\|\\d+\\|\\d+\\|", ""
         Genes.2 = str_replace_all(Genes.2,"\\d+\\|\\d+\\|\\d+\\|", ""),
         Genes.2 = str_replace_all(Genes.2,"_\\d+", "_"),
         Genes.2 = str_replace_all(Genes.2,"_.\\d+", "_"),
         Genes.2 = str_replace_all(Genes.2,"QI_", ""),
         Genes.2 = str_replace_all(Genes.2,"\\|\\d+", ""),
         Genes.2 = str_replace_all(Genes.2,"\\|--\\d+", "")
  ) %>% 
  write.table("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Gene_Duplication_Events/mazia_unique_dup_gene2.txt", col.names = F, row.names = F, quote = F, sep = '\t')
# distinct()  %>% 
# mutate(genome="EV_beadeti") %>% 

# EV_mazia_specific_OG_genes <- read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//orthofinder/Gene_Duplication_Events/mazia_unique_dup_gene2.txt", header = F, sep = '\t') %>%
#   mutate(row_names = str_extract(V1, "OG\\d+"),
#          V1 = str_remove(V1, "OG\\d+[\\s]"))
# 
# EV_mazia_specific_OG_genes %>% remove_rownames() %>% column_to_rownames(var = "row_names") %>% head()
# 
# add_rownames(EV_mazia_specific_OG_genes, var = EV_mazia_specific_OG_genes %>% 
#                mutate(row_names = str_extract(V1, "OG\\d+"),
#                       V1 = str_remove(V1, "OG\\d+[\\s]")) %>%
#                select(row_names)) 


## Read reformatted species specific OGs duplications 


library(tidyverse)

## Bedadeti

orthofinder_dups_OG_gene <- c()

# OG_Dup_Genes = c("bedadeti_unique_dup_gene1","bedadeti_unique_dup_gene2","mazia_unique_dup_gene1","mazia_unique_dup_gene2")
# OG_Dup_Genes = c("mazia_unique_dup_gene1")
OG_Dup_Genes = c("bedadeti_unique_dup_gene2","mazia_unique_dup_gene2")


getwd()

for (OG in OG_Dup_Genes){
  
  # dup_gene1 <- read.delim(paste0(OGa, ".txt"), header = F ,stringsAsFactors = FALSE) 
  unique_dup_OG <- read.delim(paste0("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Gene_Duplication_Events/", OG, "_transposed.txt"), header = F) 
  
  # bed_unique_dup_OG[1] %>% 
  #   mutate(OGs = str_extract(V1,"^OG\\d+"),
  #          OGs = str_replace_na(OGs, OGs)) %>%
  #   filter(!str_detect(V1,"^OG"), V1 !="")
  # 
  
  
  for (i in 1:ncol(unique_dup_OG)) {
    
    orthofinder_dups <- 
      unique_dup_OG[i] 
    
    
    colnames(orthofinder_dups) = c("V1")
    
    orthofinder_dups_OG_gene <-
      orthofinder_dups %>% 
      mutate(OGs = str_extract(V1,"^OG\\d+"),
             OGs = str_replace_na(OGs, OGs),
             OGs_from = OG) %>%
      filter(!str_detect(V1,"^OG"), V1 !="") %>%
      rbind(orthofinder_dups_OG_gene)
  }
  
  # orthofinder_dups_OG_gene %>% 
  #   distinct() %>% 
  #   # select(OGs) %>% distinct() %>% nrow()
  #   write.table(paste0("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//orthofinder/Gene_Duplication_Events/", OG, ".1.txt"),sep = '\t',col.names = F,row.names = F,quote = F)
}


maz_bed_dups_OG_gene <- 
  orthofinder_dups_OG_gene %>%
  mutate(OGs_from = str_remove_all(OGs_from,"_unique_dup_gene\\d"))

maz_bed_dups_OG_gene %>% 
  select(-OGs) %>%
  distinct() %>%
  group_by(OGs_from) %>%
  summarise(count_dup_OG_genes=n())
```


# UpSet 

### Plot view 
```{r upset}
library(UpSetR)
# 
EV_AED20_MAB.OG_GeneCount_Upset <-
  EV_AED20_MAB.OG_GeneCount %>%
  rename(
    
    "MA" = Musa_acuminata,
    
    "MB" = Musa_balbisiana,
    "EV (Bedadeti)" = Ensete_ventricosum_bedadeti,
    "EV (Mazia)" = Ensete_ventricosum_mazia
  ) 


EV_AED20_MAB.OG_GeneCount_Upset %>% 
  upset(
    scale.intersections="identity",
    scale.sets="identity",
    point.size=3.7,
    line.size=0.5,
    nintersect= 500,
    sets=c("EV (Bedadeti)",  "EV (Mazia)",  "MA", "MB"),
    mainbar.y.label="Number of shared\nOrthougroups",
    sets.x.label = NULL,
    order.by = c("degree"),
    text.scale = c(1.8,1.8,1.50,1.5,1.80,0.1),
    number.angles = 0,
    keep.order = TRUE)

# ggsave("../orthology/ensete_musa_v2.2/VennDiag/Upset_OGs.tiff", width=90, height=50, limitsize = FALSE)
```

### Final plot 
```{r upset_final}
EV_AED20_MAB.OG_GeneCount_Upset %>% 
  upset(
    # nsets=10,
    nintersect= 300,
    scale.intersections="identity",
    scale.sets="identity",
    point.size=3.7,
    line.size=0.5,
    # sets=c("Brassicaceae","At-pan-NLR'ome","Araport11"),
    sets=c("EV (Mazia)",  "EV (Bedadeti)",  "MA", "MB"),
    
    # "EV (Derea)",
    # sets=c("Ensete ventricosum\n(Mazia)","Ensete ventricosum\n(Bedadeti)","Ensete ventricosum\n(Derea)","Ensete ventricosum\n(Jungle Seed)","Musa acuminata","Ensete glaucum", "Musa schizocarpa","Musa balbisiana"),
    mainbar.y.label="Number of shared\nOrthougroups",
    sets.x.label = NULL,
    #      mb.ratio = c(0.65,0.35),
    order.by = c("degree"),
    # text.scale = c(1.8,1.8,1.4,1.6,1.4),
    text.scale = c(1.4,1.4,1.4,1.4,1.4,1.6),
    
    # text.scale = c("axis.y.title","axis.y.text",2.2,"col.name.side.bar","col.names",'labels.top.bar),
    # text.scale = c(2.2,2.2,1,1,2.2,2.2),
    #text.scale = c(y.title, y.axis, x.title, x.axis, sets.names, barplot.labels)
    number.angles = 0,
    # set_size.show = FALSE,
    # set_size = FALSE,
    keep.order = TRUE,
    
    
    queries = list(
      # Common OGs
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "EV (Mazia)",  "MA", "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MA", "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "EV (Mazia)",  "MA"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)","EV (Mazia)",  "MA", "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MA"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MA", "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MA", "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)", "EV (Bedadeti)",   "MA"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)", "EV (Bedadeti)",   "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)", "EV (Bedadeti)"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)", "EV (Bedadeti)" ,"MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list( "EV (Bedadeti)",   "MA"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Bedadeti)",   "MB"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Bedadeti)",   "MB"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",   "MA","MB"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)",   "MA"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)",   "MB"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)",   "MA"), 
        color = "#30ff00", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",   "MB"), 
        color = "#30ff00", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)"), 
        color = "#30ff00", 
        active = T
      ),
      
      
      
      
      # MUSA specific OGs 
      list(
        query = intersects,
        params = list("MA", "MB"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MA", "MB"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MA"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MB"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MS"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MA"), 
        color = "#d4d41c", 
        active = T
      ),
      list(
        query = intersects,
        params = list("MB"), 
        color = "#d4d41c", 
        active = T
      ),
      
      
      
      
      
      # Ensete glacucum 
      
      # list(
      #   query = intersects,
      #   params = list("MB"), 
      #   color = "#30ff00", 
      #   active = T
      # ),
      
      
      # Musa and Ensete glaucum 
      list(
        query = intersects,
        params = list( "MA", "MB"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MA"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MA","MB"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MB"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MB"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MS"),
        color = "#0BEFF1",
        active = T
      ),
      list(
        query = intersects,
        params = list( "MA"),
        color = "#0BEFF1",
        active = T
      ),
      # list(
      #   query = intersects,
      #   params = list( "MA"), 
      #   color = "#0BEFF1", 
      #   active = T
      # ),
      
      
      # Ensete OGs 
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "EV (Mazia)"), 
        color = "#155c03", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)"), 
        color = "#155c03", 
        active = T
      ),
      # list(
      #   query = intersects,
      #   params = list("EG"), 
      #   color = "#155c03", 
      #   active = T
      # ),
      list(
        query = intersects,
        params = list("EV (Mazia)"), 
        color = "#155c03", 
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)"), 
        color = "#155c03", 
        active = T
      ),
      
      # Musa and EV
      list(
        query = intersects,
        params = list("EV (Mazia)","EV (Bedadeti)",  "MA", "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MA", "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)","EV (Bedadeti)",  "MA", "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)","EV (Bedadeti)",  "MB"),
        color = "#4c01f8",
        active = T
      ),
      # list(
      #   query = intersects,
      #   params = list("EV (Mazia)","EV (Bedadeti)",  "MA"),
      #   color = "#4c01f8",
      #   active = T
      # ),
      list(
        query = intersects,
        params = list("EV (Mazia)","EV (Bedadeti)"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)","EV (Bedadeti)",  "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MA","MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MA","MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "EV (Mazia)", "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)"),
        color = "#4c01f8",
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)",  "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "EV (Mazia)", "MB"),
        color = "#4c01f8",
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Mazia)", "MB"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "MB"),
        color = "#4c01f8",
        active = T
      ),
      
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "MA"),
        color = "#4c01f8",
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)", "MB","MA"),
        color = "#4c01f8",
        active = T
      ),
      
      
      # Ensete ventrocosum specific OGs
      list(
        query = intersects,
        params = list("EV (Bedadeti)",  "EV (Mazia)"), 
        color = "#Df5286", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Bedadeti)"), 
        color = "#Df5286", 
        active = T
      ),
      list(
        query = intersects,
        params = list("EV (Mazia)"), 
        color = "#Df5286", 
        active = T
      )
    )
    # queries = list(
    #   list(
    #     query = intersects,
    #     params = list("EV (Bedadeti)"), 
    #     color = "#Df5286", 
    #     active = T
    #   )
    # )
    # 
    
  )
```







# Parse assigned OG genes 
```{r assigned OG_genes}

### final script: parse_OG_gene_list 


library(tidyverse)

## read OGs containing gene list

orthogroups_gene_list <- read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups.tsv", header = T)  

orthogroups_gene_list %>% colnames()

# read gene list count table for each OGs 
OG_GeneCount_unfiltered <- read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups.GeneCount.tsv", sep = '\t', header = T) %>% 
  dplyr::rename(
    
    "MA" = Musa_acuminata,
    "MB" = Musa_balbisiana,
    "EV_bed" = Ensete_ventricosum_bedadeti,
    "EV_maz" = Ensete_ventricosum_mazia) #%>%


# get the maximum number of genes contained in OGs for each genome 

orthogroups_gene_list_max_genes <-
  # orthogroups_gene_list %>%
  OG_GeneCount_unfiltered %>%
  # select(-Ensete_glaucum,-Ensete_ventricosum_bedadeti,-Ensete_ventricosum_mazia,
  #        -Musa_acuminata,-Musa_balbisiana,-Musa_schizocarpa) %>%
  # group_by(Orthogroup) %>%
  summarise(
    #Ensete_glaucum = max(EG),
    Ensete_ventricosum_bedadeti = max(EV_bed),
    Ensete_ventricosum_mazia = max(EV_maz),
    Musa_acuminata  = max(MA),
    Musa_balbisiana = max(MB)
    # Musa_schizocarpa = max(MS)
  ) 



Species_List = c("Musa_acuminata", "Musa_balbisiana", 
                 "Ensete_ventricosum_bedadeti", "Ensete_ventricosum_mazia")

all_assigned_OGsgenes_EV_EG_MABS <- c()


for (Each_Species in Species_List) {
  # set maximum column number for the species 
  col_range = orthogroups_gene_list_max_genes %>%
    dplyr::select(Each_Species) 
  
  
  col_range = 1:(col_range[1,1]+1)
  
  # split to genes to columns and 
  # Species_OGs <-
  all_assigned_OGsgenes_EV_EG_MABS <-
    orthogroups_gene_list %>%
    dplyr::select(Orthogroup, Each_Species) %>%
    separate(Each_Species,"seq.name", sep = "\\,")     %>%
    ## create Column IDS to combine and assign spitted genes into their respective genomes
    mutate(genome = Each_Species) %>%
    mutate(seq.name  = str_replace_na(seq.name, "") ) %>%
    dplyr::filter(seq.name !="" ) %>%
    rbind(all_assigned_OGsgenes_EV_EG_MABS)
  
  
  # Col_Names_List <- 
  #   Species_OGs %>%
  #   dplyr::select(-Orthogroup) %>% 
  #   colnames() 
  
  # for (Each_Col_Name in Col_Names_List) {
  #   
  #   all_assigned_OGsgenes_EV_EG_MABS <- 
  #     Species_OGs %>%
  #     dplyr::select(Orthogroup, Each_Col_Name) %>%
  #     dplyr::rename(pred_genes = Each_Col_Name) %>%
  #     mutate(species = Each_Species) %>%
  #     rbind(all_assigned_OGsgenes_EV_EG_MABS)
  # }
  # 
  # all_assigned_OGsgenes_EV_EG_MABS <-
  #   all_assigned_OGsgenes_EV_EG_MABS %>% 
  #   mutate(pred_genes  = str_replace_na(pred_genes, "NA") ) %>%
  #   dplyr::filter(pred_genes != "NA", pred_genes !="" ) 
  
}


all_assigned_OGsgenes_EV_EG_MABS %>% distinct(genome)


```

# Parse Unassigned OG genes 

```{r }
unassigned_prot <- "../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups_UnassignedGenes.tsv"


Orthogroups_UnassignedGenes <- read.delim(unassigned_prot, header = T, sep = '\t')

Genomes_Lists=c("Ensete_ventricosum_bedadeti",
                "Ensete_ventricosum_mazia", "Musa_acuminata", "Musa_balbisiana")

# Orthogroups_UnassignedGenes %

# Orthogroups_UnassignedGenes_EV_EG_MBS <- Orthogroups_UnassignedGenes_list

Orthogroups_UnassignedGenes_EV_EG_MBS <- c()

# for (Each_Genome in Genomes_Lists) {
#   Orthogroups_UnassignedGenes_list <- 
#     Orthogroups_UnassignedGenes %>% 
#     select(Orthogroup, Each_Genome) %>% 
#     distinct() %>%
#     # filter(Name !="") %>% 
#     # rename(seq.name=Ensete_glaucum) %>%
#     mutate(
#       # seq.name = str_extract(Name,"\\w+[^\\s]+"),
#       genome = Each_Genome) %>%
#     rename(seq.name = Each_Genome) %>%
#     rbind(Orthogroups_UnassignedGenes_list)
# }


### updated 
Orthogroups_UnassignedGenes_EV_EG_MBS  <- c()

for (Each_Genome in Genomes_Lists) {
  Orthogroups_UnassignedGenes_EV_EG_MBS <- 
    Orthogroups_UnassignedGenes %>% 
    select(Orthogroup, Each_Genome) %>% 
    distinct() %>% 
    # filter(Name !="") %>% 
    # rename(seq.name=Ensete_glaucum) %>%
    mutate(
      # seq.name = str_extract(Name,"\\w+[^\\s]+"),
      genome = Each_Genome) %>% 
    rename(seq.name = Each_Genome) %>% 
    filter(seq.name != "") %>% 
    rbind(Orthogroups_UnassignedGenes_EV_EG_MBS)
  
}


Orthogroups_UnassignedGenes_EV_EG_MBS %>% 
  filter(seq.name !="") %>% 
  group_by(genome) %>%
  summarise(count = n()) %>%
  
  ggplot(aes(x= factor(genome, levels = c("Ensete_ventricosum_bedadeti","Ensete_ventricosum_mazia", 
                                          "Ensete_glaucum","Musa_balbisiana","Musa_acuminata")),
             y=count ))+
  geom_col( position = "dodge" )+
  # scale_y_continuous(breaks =c(85,  320, 1177))+
  scale_x_discrete(breaks =c("Ensete_ventricosum_bedadeti","Ensete_ventricosum_mazia","Musa_balbisiana","Musa_acuminata"),
                   labels = c("EV (Bedadeti)","EV (Mazia)", "MB","MA"))+
  scale_y_continuous(breaks =c(1380,2398, 3311, 5955, 6418))+
  
  
  coord_flip()+
  # theme_bw() +
  labs(
    y  = " Number of unassigned orthogroups or genes",
    # x = "Genomes"
  )+
  theme(
    axis.title.x = element_markdown(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_markdown(size = 10),
    # axis.title.x
    panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 12)
    # legend.position = c(0.7,0.9)
  )

ggsave("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB//Orthogroups/VennDiag/Unassigned_OG_genes.tiff", width=7, height=7)
```





# Orthogroups by type: Multi-single-unique-unassinged 
## prepare data
```{r ortho_type}
unfiltered_EV_AED20_MAB.OG_GeneCount <-
  read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups.GeneCount.tsv", sep = '\t', header = T) %>%
  bind_rows(
    read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>% 
      mutate(
        
        Ensete_ventricosum_mazia = case_when(str_detect(Ensete_ventricosum_mazia,"\\w+") ~ 1,
                                             TRUE ~ 0),
        Ensete_ventricosum_bedadeti = case_when(str_detect(Ensete_ventricosum_bedadeti, "\\w+") ~ 1,
                                                TRUE ~ 0),
        Musa_acuminata = case_when(str_detect(Musa_acuminata, "\\w+") ~ 1,
                                   TRUE ~ 0),
        Musa_balbisiana = case_when(str_detect(Musa_balbisiana, "\\w+") ~ 1,
                                    TRUE ~ 0)
      ))


single_copy_ortholgos <- 
  unfiltered_EV_AED20_MAB.OG_GeneCount %>%
  filter(
    
    Musa_acuminata == 1,
    Musa_balbisiana == 1,
    Ensete_ventricosum_bedadeti == 1,
    Ensete_ventricosum_mazia ==1
  ) %>%
  mutate (
    single_copy_OGs = "yes"
  )

# summarise(
#   Ensete_glaucum = sum(Ensete_glaucum), 
#   Musa_acuminata = sum(Musa_acuminata),
#   Musa_balbisiana = sum(Musa_balbisiana),
#   Musa_schizocarpa = sum(Musa_schizocarpa),
#   Ensete_ventricosum_bedadeti = sum(Ensete_ventricosum_bedadeti),
#   Ensete_ventricosum_mazia = sum(Ensete_ventricosum_mazia)
# ) %>% 
# pivot_longer( cols = Ensete_glaucum:Ensete_ventricosum_mazia,
#               values_to = "value")%>%
# mutate(Orthogroup = "Single-copy orthologs")

unassigned_OGs <-
  unfiltered_EV_AED20_MAB.OG_GeneCount %>%
  left_join(
    single_copy_ortholgos
  ) %>% 
  mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
  filter(single_copy_OGs == "NA") %>%
  filter(
    
    Musa_acuminata == 1,
    Musa_balbisiana == 0,
    Ensete_ventricosum_bedadeti == 0,
    Ensete_ventricosum_mazia ==0
  ) %>% 
  rbind(
    
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana == 1,
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia ==0
      ),
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana == 0,
        Ensete_ventricosum_bedadeti == 1,
        Ensete_ventricosum_mazia ==0
      ),
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana == 0,
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia ==1
      )
  ) %>% 
  mutate(unassigned_OG = "yes") %>%
  select(-single_copy_OGs)


single_copy_ortholgos %>% head()


unique_OGs <-
  unfiltered_EV_AED20_MAB.OG_GeneCount %>%
  left_join(
    single_copy_ortholgos
  ) %>% 
  mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
  filter(single_copy_OGs == "NA") %>%
  filter(
    
    Musa_acuminata > 1,
    Musa_balbisiana == 0,
    Ensete_ventricosum_bedadeti == 0,
    Ensete_ventricosum_mazia ==0
  ) %>% 
  rbind(
    
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana > 1,
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia ==0
      ),
    
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana == 0,
        Ensete_ventricosum_bedadeti > 1,
        Ensete_ventricosum_mazia ==0
      ),
    unfiltered_EV_AED20_MAB.OG_GeneCount %>%
      left_join(
        single_copy_ortholgos
      ) %>% 
      mutate(single_copy_OGs = str_replace_na(single_copy_OGs,"NA")) %>%
      filter(single_copy_OGs == "NA") %>%
      filter(
        
        Musa_acuminata == 0,
        Musa_balbisiana == 0,
        Ensete_ventricosum_bedadeti == 0,
        Ensete_ventricosum_mazia > 1
      )
  ) %>% 
  mutate(unique_OGs = "yes") %>%
  select(-single_copy_OGs) 

## multi-copy OGs 

multi_copy_OGs <-
  unfiltered_EV_AED20_MAB.OG_GeneCount %>% 
  left_join(
    single_copy_ortholgos
  ) %>%
  left_join(
    unique_OGs
  ) %>% 
  left_join(
    unassigned_OGs
  ) %>% 
  mutate(
    single_copy_OGs = str_replace_na(single_copy_OGs,"NA"),
    unique_OGs = str_replace_na(unique_OGs,"NA"),
    unassigned_OG = str_replace_na(unassigned_OG,"NA")
  ) %>% 
  filter(
    single_copy_OGs == "NA" ,
    unique_OGs == 'NA' ,
    unassigned_OG == "NA"
  ) %>%
  distinct() %>%
  select(-single_copy_OGs, -unique_OGs, -unassigned_OG, -Total)

# multi_copy_OGs %>% head()

# unique_OGs
# single_copy_ortholgos
# unassigned_OGs

multi_single_unique_unassigned_OG_summary <-
  multi_copy_OGs %>%
  summarise(
    # Ensete_glaucum = sum(Ensete_glaucum),
    Musa_acuminata = sum(Musa_acuminata),
    Musa_balbisiana = sum(Musa_balbisiana),
    # Musa_schizocarpa = sum(Musa_schizocarpa),
    Ensete_ventricosum_bedadeti = sum(Ensete_ventricosum_bedadeti),
    Ensete_ventricosum_mazia = sum(Ensete_ventricosum_mazia)
  ) %>% 
  pivot_longer( cols = Musa_acuminata:Ensete_ventricosum_mazia,
                values_to = "value")%>%
  mutate(Orthogroup = "Multi-copy orthologs") %>%
  rbind(
    single_copy_ortholgos %>%
      summarise(
        # Ensete_glaucum = sum(Ensete_glaucum),
        Musa_acuminata = sum(Musa_acuminata),
        Musa_balbisiana = sum(Musa_balbisiana),
        # Musa_schizocarpa = sum(Musa_schizocarpa),
        Ensete_ventricosum_bedadeti = sum(Ensete_ventricosum_bedadeti),
        Ensete_ventricosum_mazia = sum(Ensete_ventricosum_mazia)
      ) %>%
      pivot_longer( cols = Musa_acuminata:Ensete_ventricosum_mazia,
                    values_to = "value")%>%
      mutate(Orthogroup = "Single-copy orthologs"),
    unique_OGs %>%
      summarise(
        # Ensete_glaucum = sum(Ensete_glaucum),
        Musa_acuminata = sum(Musa_acuminata),
        Musa_balbisiana = sum(Musa_balbisiana),
        # Musa_schizocarpa = sum(Musa_schizocarpa),
        Ensete_ventricosum_bedadeti = sum(Ensete_ventricosum_bedadeti),
        Ensete_ventricosum_mazia = sum(Ensete_ventricosum_mazia)
      ) %>%
      pivot_longer( cols = Musa_acuminata:Ensete_ventricosum_mazia,
                    values_to = "value")%>%
      mutate(Orthogroup = "Unique"),
    unassigned_OGs %>%
      summarise(
        # Ensete_glaucum = sum(Ensete_glaucum),
        Musa_acuminata = sum(Musa_acuminata),
        Musa_balbisiana = sum(Musa_balbisiana),
        # Musa_schizocarpa = sum(Musa_schizocarpa),
        Ensete_ventricosum_bedadeti = sum(Ensete_ventricosum_bedadeti),
        Ensete_ventricosum_mazia = sum(Ensete_ventricosum_mazia)
      ) %>%
      pivot_longer( cols = Musa_acuminata:Ensete_ventricosum_mazia,
                    values_to = "value")%>%
      mutate(Orthogroup = "Unassigned")
  ) 
```


## Assigned_unassigned_combined
```{r all OGs}

assigned_unassigned_OGs_EV_EG_MABS <-
  multi_copy_OGs %>% 
  mutate (
    Ensete_ventricosum_bedadeti = case_when(Ensete_ventricosum_bedadeti > 1 ~ "bedadeti" ),
    Ensete_ventricosum_mazia= case_when(Ensete_ventricosum_mazia > 1 ~ "mazia"),
    Musa_acuminata= case_when(Musa_acuminata > 1 ~ "musa_ac"),
    Musa_balbisiana= case_when(Musa_balbisiana > 1 ~ "musa_ba")) %>%
  bind_rows(
    unique_OGs  %>%
      mutate (
        Ensete_ventricosum_bedadeti = case_when(Ensete_ventricosum_bedadeti >=1 ~ "bedadeti" ),
        Ensete_ventricosum_mazia= case_when(Ensete_ventricosum_mazia >=1 ~ "mazia"),
        Musa_acuminata= case_when(Musa_acuminata >=1 ~ "musa_ac"),
        Musa_balbisiana= case_when(Musa_balbisiana >=1 ~ "musa_ba")), 
    unassigned_OGs %>%
      mutate (
        Ensete_ventricosum_bedadeti = case_when(Ensete_ventricosum_bedadeti ==1 ~ "bedadeti" ),
        Ensete_ventricosum_mazia= case_when(Ensete_ventricosum_mazia ==1 ~ "mazia"),
        Musa_acuminata= case_when(Musa_acuminata ==1 ~ "musa_ac"),
        Musa_balbisiana= case_when(Musa_balbisiana ==1 ~ "musa_ba")),
    single_copy_ortholgos %>%
      mutate (
        Ensete_ventricosum_bedadeti = case_when(Ensete_ventricosum_bedadeti ==1 ~ "bedadeti" ),
        Ensete_ventricosum_mazia= case_when(Ensete_ventricosum_mazia ==1 ~ "mazia"),
        Musa_acuminata= case_when(Musa_acuminata ==1 ~ "musa_ac"),
        Musa_balbisiana= case_when(Musa_balbisiana ==1 ~ "musa_ba"))) 

assigned_unassigned_OGs_EV_EG_MABS_v1 <- 
  assigned_unassigned_OGs_EV_EG_MABS %>% 
  pivot_longer( cols = Musa_acuminata:Ensete_ventricosum_mazia,
                values_to = "value") %>%
  select(Orthogroup,name,value) %>%
  mutate(value = str_replace_na(value, "NA")) %>%
  filter(value != "NA") %>%
  distinct()

# assigned_unassigned_OGs_EV_EG_MABS_v1 %>%
#   rename(genome = name) %>%
#   left_join(
#     
#     Orthogroups_UnassignedGenes_EV_EG_MBS %>% 
#       head() %>%
#       bind_rows(
#         all_assigned_OGsgenes_EV_EG_MABS) %>% 
#       # rename(#seq.name = pred_genes, 
#       #        genome = species)) %>%
#       distinct()) %>% tail()


EV_AED20_MAB.OG_GeneCount %>%
  filter(Orthogroup == "OG0016389")

assigned_unassigned_OGs_genes_EV_EG_MAB <- 
  read.delim("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/Orthogroups/Orthogroups.txt", sep = ' ', header = F) %>%  
  pivot_longer( cols = V2:V1751,
                values_to = "seq.name") %>%
  mutate(V1 = str_remove(V1,":"),
         genome = case_when (str_detect(seq.name,"EVBD") ~ "bedadeti",
                             str_detect(seq.name,"EVMZ") ~ "mazia",
                             # str_detect(seq.name, "Eg\\d+") ~ "ensete_gl",
                             str_detect(seq.name, "Mba") ~ "musa_ba",
                             str_detect(seq.name, "Macma") ~ "musa_ac",
                             # str_detect(seq.name) ~ "musa_sc"
         )) %>%
  rename(Orthogroup = V1) %>% 
  select(Orthogroup, seq.name,genome) %>%
  mutate(genome = str_replace_na(genome, "NA")) %>%
  filter(genome != "NA") %>%
  
  # bind unassigned OGs
  
  bind_rows(Orthogroups_UnassignedGenes_EV_EG_MBS)



assigned_unassigned_OGs_genes_EV_EG_MAB %>%
  tail()
```



## plot

```{r plot}
multi_single_unique_unassigned_OG_summary.v1 <- 
  multi_single_unique_unassigned_OG_summary %>%
  # filter(name == "Ensete_ventricosum_mazia")
  mutate(name = str_replace_all(name,c("Ensete_ventricosum_bedadeti" = "*E. ventricosum* (Bedadeti)",
                                       "Ensete_ventricosum_mazia" = "*E. ventricosum* (Mazia)",
                                       "Musa_balbisiana" = "*M. balbisiana*",
                                       "Musa_acuminata" = "*M. acuminata*"))) 



multi_single_unique_unassigned_OG_summary.v1$name <- factor(multi_single_unique_unassigned_OG_summary.v1$name, 
                                                            levels = c("*E. ventricosum* (Mazia)", 
                                                                       "*E. ventricosum* (Bedadeti)",
                                                                       "*M. balbisiana*",  
                                                                       "*M. acuminata*"))

multi_single_unique_unassigned_OG_summary.v1 %>%
  
  ggplot(aes(name,value, fill = Orthogroup)) +
  geom_col()+
  # scale_y_continuous(breaks =c(0, 2239, 12890, 26000,39000))+
  # scale_x_discrete(breaks =c( "MB", "MA", "EV_bedadeti", "EV_mazia"),
  #                  labels = c("MB","MA",  "EV (Bedadeti)","EV (Mazia)"))+
  
  coord_flip()+
  # theme_bw() +
  labs(
    y =" Number of genes",
    # x = "Genomes"
  )+
  theme(
    axis.title.x = element_markdown(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_markdown(size = 12),
    # axis.title.x
    panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_blank(),
    # legend.text = element_markdown(size = 12),
    legend.text = element_markdown(size=10.8),
    legend.key.size = unit(0.2,"cm"), 
    legend.position = "top"
  )

ggsave("../../reannotation_analysis/mazia_reannotation/orthofinder/EV_AED20_MAB/OG_EV_MAB_summary.tiff", width=7, height=4)

```

```{r save}
save.image(file = "../EV_MAB/OrthoFinder.EV.MAB.RData")
```