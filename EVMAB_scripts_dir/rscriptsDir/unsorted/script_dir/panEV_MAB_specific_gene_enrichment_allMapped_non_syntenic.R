## Updated data
EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  mutate(syn_cov = str_replace_na(syn_cov,"NA"),
         synteny = str_replace_na(synteny,"NA")) %>%
  # filter(synteny=="NA") %>%
  filter(gene_id=="EVBD.1.000038") 

  distinct(syn_cov)
  head()

## previoiuse with reads mapping only dasta
# EV_MAB_reads_coverage_geneID_GO_terms %>%
  # head()

## Vendiagram plot for EV PAV genes 
### For Mazia PAV genes 

EV_MAB_reads_coverage_geneID_GO_terms_wide <-
  EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  # filter(str_detect(GO_terms,"GO:")) %>%
  # head()
  # dplyr::select(gene_id, fraction_overlap) %>%
  mutate(reads_genome = str_replace(reads_genome,"EV genotypes", "EV"),
         reads_genome = str_replace(reads_genome,"A-sub genome banana genotypes", "MA"),
         reads_genome = str_replace(reads_genome,"B-sub genome banana genotypes", "MB"),
         
         # For non-syntenic genes set the coverage value to 1 meaning that they showed some sort of nucmer or lastz alignment based synteny  
         # cov = str_replace_na(cov,"1")
         # cov = as.numeric(cov)
         ) %>% 
    # head()
  # distinct(reads_used)
  dplyr::filter(reads_used != "qfiltered") %>%
  dplyr::filter(gene_repeat == "gene")
  # head()
EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  mutate(syn_cov = str_replace_na(syn_cov,"NA"),
         synteny = str_replace_na(synteny,"NA")) %>%
  filter(syn_cov=="NA") %>%
  filter(synteny == "NA") %>%
head()
  


EV_MAB_reads_coverage_geneID_GO_terms_wide <-
  EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
    mutate(syn_cov = str_replace_na(syn_cov,"NA"),
           synteny = str_replace_na(synteny,"NA")) %>%
  dplyr::select(gene_id, ref_genome, query, reads_genome,synteny,GO_terms,fraction_overlap,syn_cov) %>%
    # filter(synteny == "syntenic") %>%
    # filter(syn_cov < 0.25) %>%
    mutate(
      # set fraction_coverage to 0 if the is non-syntenic and the initial reads coverage is less than 0.25%
      syn_frac_cov = case_when(
        
        # This could be de novo evolved genes who are not syntenic and do not have reads support
        synteny == "non_syntenic" & fraction_overlap < 0.25 ~ "0",
        
        # This could be highly diverged genes
        synteny == "NA" & fraction_overlap < 0.25 ~ "0",
        
        # synteny == "NA" & fraction_overlap < 0.25 & syn_cov < 0.25 ~ "0",
        
        # Genes that show syntenic relationship but with lower syntenic and reads mapping coverage.
        # This may be a small proteins or a small conserved protein domain commmon in Musaceae 
        # (May be not and it's not the interest of the analyis we are looking at)  
        synteny == "syntenic" & fraction_overlap < 0.25 & syn_cov < 0.25 ~ "0",
        
      
        # use fraction_overlap for intera-species reads mapping
        synteny == "self_reads" & fraction_overlap < 0.25  ~ "0",
        
        TRUE ~ "1"),
      # test = case_when(str_detect(synteny,"non_syntenic") ~ "0",
      #                              TRUE ~ fraction_overlap)
      syn_frac_cov = as.numeric(syn_frac_cov)
    ) %>% 
    select(-syn_cov,-fraction_overlap,-synteny ) %>%
    group_by(gene_id,ref_genome, reads_genome,GO_terms) %>%
    summarise(syn_frac_cov = max(syn_frac_cov)) %>%
    ungroup() %>%
    pivot_wider( names_from = reads_genome, values_from = syn_frac_cov) %>% 
    as.data.frame ()


EV_MAB_reads_coverage_geneID_GO_terms_wide %>% 
  filter(ref_genome=="ev_bedadeti") %>%
  # mutate(EV = str_replace_na(EV, "NA")) %>%
  # filter(EV == 0) %>%
  filter(gene_id=="EVBD.1.000038") 

  head()


EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  mutate(syn_cov = str_replace_na(syn_cov,"NA"),
         synteny = str_replace_na(synteny,"NA")) %>%
  filter(ref_genome=="ev_bedadeti") %>%
  filter(gene_id=="EVBD.1.000038") 
  
EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
  filter(end=="125076") %>% 
  as.data.frame()

## For Bedadeti genome 

EV_MAB_reads_coverage_geneID_GO_terms_wide %>% 
  as.data.frame() %>% 
  filter(ref_genome=="ev_bedadeti") %>% 
  # filter(query == "panMAB_reads") %>% 
  filter(MA ==1) %>% 
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_bedadeeti_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="ev_bedadeti") %>%
  # filter(query == "panMAB_reads") %>%
  filter(MB == 1) %>%
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_bedadeeti_panMB.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="ev_bedadeti") %>%
  # filter(query == "panEV_reads") %>%
  filter(EV == 1) %>%
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_bedadeeti_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)


# Figure of main venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

# OG_list_VennDiag_EV_MAB[["EV (Bedadeti)"]] <- read.table("../monocot_blastp/advanced_orthofinder/OrthogroupsVennDiag/Ensete_ventricosum_bedadeti.OG.txt", header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(DATA_DIR, "EV_bedadeeti_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(DATA_DIR, "EV_bedadeeti_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(DATA_DIR, "EV_bedadeeti_panMB.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

# Plot Ven diagra 
library(VennDiagram)
library(ggVennDiagram)

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

EV_bedadeti_vs_MAB_vendiag <-
  ggVennDiagram(x = Gene_list_VennDiag_PanEV_MAB, 
                edge_size = 1.5,
                size = 6,
                set_size = 0,
                label_size = 4.3,
                set_color = c("#339933","#FF7600", "#000099"))+
  # scale_fill_gradient(low="white", high = "red", limits = c(0, 40000)) +
  scale_fill_gradient(low="white", high = "white", limits = c(0, 40000)) +
  
  labs(fill="Count",
       title =NULL
       # "Presence / absence of EV (Mazia) predicted genes\n in EV, MA and MB Pan-reads"
  ) +
  # ggtitle("Shared Orthogroups (OGs)") +
  theme(#plot.title = element_text(hjust = 0.6),
    legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .4))+
  scale_y_continuous(expand = expansion(mult = .09))
  
  ## looking into EV genes (83) that lack reads coverage from enset genotypes 
  
  # EV genes that do not have EV reads support  
  venn_data_process$regionData$item[[6]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    head() %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    distinct()

  ## Collect gene sets 
  
  ### Bedadeti specific 

  EV_bedadeti_specific <- 
  venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  
  EV_bedadeti_MAB_core <-
    venn_data_process$regionData$item[[7]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
      # distinct(gene_id) %>%
      # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  
  EV_bedadeti_MA_specific <-
    venn_data_process$regionData$item[[4]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
    
  
  EV_bedadeti_MB_specific <-
    venn_data_process$regionData$item[[5]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  
  EV_bedadeti_noEVReadSupport_Butsyntenic2MAB <-
    venn_data_process$regionData$item[[6]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  

  
  
  
## For mazia genome 
EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="ev_mazia") %>%
  filter(query == "panMAB_reads") %>%
  filter(MA ==1) %>%
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_mazia_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="ev_mazia") %>%
  filter(query == "panMAB_reads") %>%
  filter(MB == 1) %>%
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_mazia_panMB.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="ev_mazia") %>%
  filter(query == "panEV_reads") %>%
  filter(EV == 1) %>%
  select(gene_id) %>%
  distinct() %>% 
  write.table(paste0(DATA_DIR, "EV_mazia_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome == "ev_mazia") %>%
  distinct(gene_id) %>%
  distinct() %>%
  nrow()


# Figure of main venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

# OG_list_VennDiag_EV_MAB[["EV (Bedadeti)"]] <- read.table("../monocot_blastp/advanced_orthofinder/OrthogroupsVennDiag/Ensete_ventricosum_bedadeti.OG.txt", header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(DATA_DIR, "EV_mazia_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(DATA_DIR, "EV_mazia_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(DATA_DIR, "EV_mazia_panMB.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

# Plot Ven diagra 
# library(VennDiagram)
# library(ggVennDiagram)

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

EV_mazia_vs_MAB_vendiag <-
  ggVennDiagram(x = Gene_list_VennDiag_PanEV_MAB, 
                edge_size = 1.5,
                size = 6,
                set_size = 0,
                label_size = 4.3,
                set_color = c("#339933","#FF7600", "#000099"))+
  # coord_flip()+
  # scale_fill_gradient(low="white", high = "red", limits = c(0, 40000)) +
  scale_fill_gradient(low="white", high = "white", limits = c(0, 40000)) +
  
  labs(fill="Count",
       title =NULL
       # "Presence / absence of EV (Mazia) predicted genes\n in EV, MA and MB Pan-reads"
  ) +
  # ggtitle("Shared Orthogroups (OGs)") +
  theme(#plot.title = element_text(hjust = 0.6),
    legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .4))+
  scale_y_continuous(expand = expansion(mult = .09))


## Collect gene sets 

EV_mazia_specific <-
  venn_data_process$regionData$item[[1]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

EV_bedadeti_MAB_core
EV_mazia_MAB_core <-
  venn_data_process$regionData$item[[7]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(GO_terms,gene_id)  %>%
  distinct() 


EV_mazia_MA_specific <-
  venn_data_process$regionData$item[[4]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(GO_terms,gene_id)  %>%
  distinct() 


EV_mazia_MB_specific <-
  venn_data_process$regionData$item[[5]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(GO_terms,gene_id)  %>%
  distinct() 


EV_mazia_noEVReadSupport_Butsyntenic2MAB <-
  venn_data_process$regionData$item[[6]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(GO_terms,gene_id)  %>%
  distinct() 


## For Musa acuminata
# EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
#   head()

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="musa_acuminata") %>%
  filter(query == "panEV_reads") %>%
  filter(EV ==1) %>%
  select(gene_id) %>%
  write.table(paste0(DATA_DIR, "MA_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="musa_acuminata") %>%
  filter(query == "panMA_reads") %>%
  filter(MA == 1) %>%
  select(gene_id) %>%
  write.table(paste0(DATA_DIR, "MA_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)


# Figure of main venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(DATA_DIR, "MA_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(DATA_DIR, "MA_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

MA_vs_EV_vendiag <-
  ggVennDiagram(x = Gene_list_VennDiag_PanEV_MAB, 
                edge_size = 1.5,
                size = 6,
                set_size = 0,
                label_size = 4.3,
                set_color = c("#339933","#FF7600"))+
  coord_flip() +
  
  # scale_fill_gradient(low="white", high = "red", limits = c(0, 40000)) +
  scale_fill_gradient(low="white", high = "white", limits = c(0, 40000)) +
  
  labs(fill="Count",
       title =NULL
       # "Presence / absence of EV (Mazia) predicted genes\n in EV, MA and MB Pan-reads"
  ) +
  # ggtitle("Shared Orthogroups (OGs)") +
  theme(#plot.title = element_text(hjust = 0.6),
    legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .4))+
  scale_y_continuous(expand = expansion(mult = .09))


  ## looking into MA genes (26) that lack reads coverage from MA genotypes 
  
  # MA genes that do not have MA reads support  
  venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    head() %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    distinct()
  
  
  
  ## collect gene sets
  MA_specific <-
    venn_data_process$regionData$item[[2]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  EV_bedadeti_MAB_core
  EV_mazia_MAB_core
  MA_EV_core <-
    venn_data_process$regionData$item[[3]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  MA_noMAReadSupport_Butsyntenic2EV <-
    venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  
  
## For Musa balbisiana
EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  distinct(ref_genome)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="musa_balbisiana") %>%
  filter(query == "panEV_reads") %>%
  filter(EV ==1) %>%
  select(gene_id) %>%
  write.table(paste0(DATA_DIR, "MB_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  filter(ref_genome=="musa_balbisiana") %>%
  filter(query == "panMB_reads") %>%
  filter(MB == 1) %>%
  select(gene_id) %>%
  write.table(paste0(DATA_DIR, "MB_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)

# Figure of main venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(DATA_DIR, "MB_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(DATA_DIR, "MB_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

MB_vs_EV_vendiag <-
  ggVennDiagram(x = Gene_list_VennDiag_PanEV_MAB, 
                edge_size = 1.5,
                size = 6,
                set_size = 0,
                label_size = 4.3,
                set_color = c("#339933", "#000099"))+
  coord_flip() +
  
  scale_fill_gradient(low="white", high = "white", limits = c(0, 40000)) +
  labs(fill="Count",
       title =NULL
  ) +
  theme(#plot.title = element_text(hjust = 0.6),
    legend.position = "none")+
  scale_x_continuous(expand = expansion(mult = .4))+
  scale_y_continuous(expand = expansion(mult = .09))
  
  
  ## looking into MB genes (26) that lack reads coverage from MB genotypes 
  
  # MB genes that do not have MN reads support  
  venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    head() %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    distinct() 
  
  MB_specific <-
    venn_data_process$regionData$item[[2]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  

  MB_EV_core <-
    venn_data_process$regionData$item[[3]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  MB_noMBReadSupport_Butsyntenic2EV <-
    venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_wide
    ) %>%
    # distinct(gene_id) %>%
    # nrow()
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  
  
  

library(patchwork)

EV_mazia_vs_MAB_vendiag+EV_bedadeti_vs_MAB_vendiag
ggsave(paste0(DATA_DIR, "plot_out/EV_mazia_bedadeti_vendiag.allMapped.synteny.pdf"), width=10, height=4)

MA_vs_EV_vendiag+MB_vs_EV_vendiag
ggsave(paste0(DATA_DIR, "plot_out/MA_MB_vendiag.allMapped.synteny.pdf"), width=6, height=4)


(EV_mazia_vs_MAB_vendiag/MA_vs_EV_vendiag)
# /(MA_vs_EV_vendiag|MB_vs_EV_vendiag)

ggsave(paste0(DATA_DIR, "plot_out/EV_mazia_top_vs_MA_bottom_genes_vendiag.allMapped.synteny.pdf"), width=9, height=9)


(EV_bedadeti_vs_MAB_vendiag/MB_vs_EV_vendiag)
# /(MA_vs_EV_vendiag|MB_vs_EV_vendiag)

ggsave(paste0(DATA_DIR, "plot_out/EV_bedadeti_top_vs_MB_bottom_genes_vendiag.allMapped.synteny.pdf"), width=9, height=9)

save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.RData"))
save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.02.RData"))

load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.RData"))
