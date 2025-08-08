# import the final species specifc genes list 
result_dir="Z:/sadik/evmab_dir/data_dir/EV_MAB_compgeno/EV_MAB_compgeno/result_dir/"
SPECIES_SPECIFIC_CDS_DIR=paste0(DATA_DIR, "synteny_dir/final_species_specific_genes/")
list.files(path = SPECIES_SPECIFIC_CDS_DIR, pattern = ".bed")
list.files(path = result_dir, pattern = "cov.bed")

# EV mazia CDS >25% panEV CDS coverage 

EV_maiza_panEV_reads_specific <-
EV_MAB_reads_coverage %>% 
  filter(reads_used != "qfiltered" ) %>%
  filter(reads_genome != "S-sub genome banana genotypes",
         reads_genome != "T-sub genome banana genotypes") %>%
  # distinct(query) 
  
  mutate(fraction_overlap = round(fraction_overlap,2)) %>%
  filter(ref_genome == "ev_mazia") %>%
  # distinct(query)
  filter(query == "panEV_reads") %>%
  select(seq.name,start,end,fraction_overlap) %>%
  group_by(seq.name,start,end) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>% 
  filter(fraction_overlap >= 0.25) 

# EV_maiza_CDS_lacking_panEV_reads_support <-
  EV_MAB_reads_coverage %>% 
  filter(reads_used != "qfiltered" ) %>%
  filter(reads_genome != "S-sub genome banana genotypes",
         reads_genome != "T-sub genome banana genotypes") %>%
  # distinct(query) 
  
  mutate(fraction_overlap = round(fraction_overlap,2)) %>%
  filter(ref_genome == "ev_mazia") %>%
  # distinct(query)
  filter(query == "panEV_reads") %>%
  select(seq.name,start,end,fraction_overlap) %>%
  group_by(seq.name,start,end) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>% 
  filter(fraction_overlap < 0.25) %>%
    write.table(paste0(result_dir,"EV_mazia_CDS_lacking_panEV_reads_support.bed"),
                col.names = F,
                row.names = F,
                quote = F,
                sep = '\t')

  # Create ouptut directory 
  dir.create(paste0(result_dir,"venndiags"))
  venndiagOut <- paste0(result_dir,"venndiags/")

  # For Bedadeti genome against MA reads mapping and/or CDS alignment
  
  ## Against MA reads mapping and/or CDS alignment
  read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"bedadeti_specific_against_panma_nucmer_lastz.MA.CDS.bed"),
             header = F) %>%
              rename(
                seqid = V1,
                start = V2,
                end = V3,
                lastzCov = V4
              ) %>%
    distinct() %>%
    left_join(
      # Join with list of bedadeti specific genes from panma reads
      read.delim(paste0(result_dir,"filtered_ev_bedadeti_cds_less_25_per_panma_cov.bed"), header = F) %>%
        mutate(source = "panMA_reads") %>%
        rename(seqid = V1,
               start = V2,
               end = V3)) %>%
    mutate(source = str_replace_na(source, "NA")) %>%
    # join with CDS gene IDs and gene ontology (GO) terms. 
    full_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_bedadeti") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov < 24.4) %>%
    distinct(gene_id,lastzCov) %>%
    # Join all gene ids of ref genome 
    full_join(EV_MAB_reads_coverage_geneID_GO_terms %>% 
    filter(ref_genome == "ev_bedadeti") %>% 
    distinct(gene_id, ref_genome) %>%
    mutate(gene_id_all = gene_id)) %>% 
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # take genes that were not included in the lastz filtered alignment 
    filter(lastzCov == "NA") %>%
    # tail()
    distinct(gene_id) %>%
    write.table(paste0(venndiagOut, "EV_bedadeeti_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  # Against MB reads mapping and/or CDS alignment
  read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"bedadeti_specific_against_panmb_nucmer_lastz.MB.CDS.bed"),
             header = F) %>%
    rename(
      seqid = V1,
      start = V2,
      end = V3,
      lastzCov = V4
    ) %>%
    distinct() %>% 
    left_join(
    
    # Join with list of bedadeti specific genes from panmb reads
    read.delim(paste0(result_dir,"filtered_ev_bedadeti_cds_less_25_per_panmb_cov.bed"), header = F) %>%
      mutate(source = "panMB_reads") %>%
      rename(seqid = V1,
             start = V2,
             end = V3)) %>%
    mutate(source = str_replace_na(source, "NA")) %>%
    # join with CDS gene IDs and gene ontology (GO) terms. 
    full_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_bedadeti") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov < 24.4 ) %>%
    distinct(gene_id,lastzCov) %>%
    # distinct(gene_id) %>%
    # nrow()
    # Join all gene ids of ref genome 
    full_join(EV_MAB_reads_coverage_geneID_GO_terms %>% 
                filter(ref_genome == "ev_bedadeti") %>% 
                distinct(gene_id, ref_genome) %>%
                mutate(gene_id_all = gene_id)) %>% 
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # take genes that were not included in the lastz filtered alignment 
    filter(lastzCov == "NA") %>%
    # tail()
    distinct(gene_id) %>%
    # nrow()
    write.table(paste0(venndiagOut, "EV_bedadeeti_panMB.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  # Against EV reads mapping 
  EV_MAB_reads_coverage %>%
    filter(reads_used != "qfiltered" ) %>%
    filter(reads_genome != "S-sub genome banana genotypes",
           reads_genome != "T-sub genome banana genotypes") %>%
    # distinct(query) 
    
    mutate(fraction_overlap = round(fraction_overlap,2)) %>%
    filter(ref_genome == "ev_bedadeti") %>%
    # distinct(query)
    filter(query == "panEV_reads") %>%
    select(seq.name,start,end,fraction_overlap) %>%
    group_by(seq.name,start,end) %>%
    summarise(fraction_overlap = max(fraction_overlap)) %>%
    ungroup() %>% 
    filter(fraction_overlap >= 0.25) %>%
    distinct(seq.name, start, end) %>%
    rename(
      seqid = seq.name
    ) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_bedadeti") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    distinct(gene_id) %>%
    write.table(paste0(venndiagOut, "EV_bedadeeti_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  

  ## Generate venn diagrams
  
  Gene_list_VennDiag_PanEV_MAB <- list()
  
  # OG_list_VennDiag_EV_MAB[["EV (Bedadeti)"]] <- read.table("../monocot_blastp/advanced_orthofinder/OrthogroupsVennDiag/Ensete_ventricosum_bedadeti.OG.txt", header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(venndiagOut, "EV_bedadeeti_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(venndiagOut, "EV_bedadeeti_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(venndiagOut, "EV_bedadeeti_panMB.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  
  # Plot Ven diagra 
  library(VennDiagram)
  library(ggVennDiagram)
  
  venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
  venn_data_process <- process_data(venn_data)
  
  EV_bedadeti_vs_MAB_vendiag.v2 <-
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
venn_data_process$regionData$item[[7]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  # head() %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
    # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
    #   select(-GO_terms) %>%
    #   filter(reads_used != "qfiltered")  %>%
    #   filter(gene_repeat == "gene")
  ) %>%
  distinct(gene_id) %>%
  nrow ()


## Collect gene sets 

EV_bedadeti_specific.v2 <-
  venn_data_process$regionData$item[[1]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  select(GO_terms,gene_id)  %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

EV_bedadeti_MAB_core.v2 <-
  venn_data_process$regionData$item[[7]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

EV_bedadeti_MA_specific.v2 <-
  venn_data_process$regionData$item[[4]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

EV_bedadeti_MB_specific.v2 <-
  venn_data_process$regionData$item[[5]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

# EV_bedadeti_noEVReadSupport_Butsyntenic2MAB <-
  venn_data_process$regionData$item[[6]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 


  # For EV (Mazia) genome against MA reads mapping and/or CDS alignment
  
  ## Against MA reads mapping and/or CDS alignment
  read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"mazia_specific_against_panma_nucmer_lastz.MA.CDS.bed"),
             header = F) %>%
    rename(
      seqid = V1,
      start = V2,
      end = V3,
      lastzCov = V4
    ) %>%
    distinct() %>%
    # # nrow()
    # # join with CDS gene IDs and gene ontology (GO) terms. 
    # full_join(
    #   EV_MAB_reads_coverage_geneID_GO_terms %>% 
    #     filter(ref_genome == "ev_mazia") %>% 
    #     select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
    #     distinct() 
    # ) %>%
    # mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # filter(lastzCov == "NA") %>%
    # distinct(gene_id) %>%
    left_join(
      # Join with list of mazia specific genes from panma reads
      read.delim(paste0(result_dir,"filtered_ev_mazia_cds_less_25_per_panma_cov.bed"), header = F) %>%
        mutate(source = "panMA_reads") %>%
        rename(seqid = V1,
               start = V2,
               end = V3)) %>%
    mutate(source = str_replace_na(source, "NA")) %>%
    # join with CDS gene IDs and gene ontology (GO) terms. 
    full_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_mazia") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov < 24.4) %>%
    distinct(gene_id,lastzCov) %>%
    # distinct(gene_id) %>%
    # nrow()
    # Join all gene ids of ref genome 
    full_join(EV_MAB_reads_coverage_geneID_GO_terms %>% 
                filter(ref_genome == "ev_mazia") %>% 
                distinct(gene_id, ref_genome) %>%
                mutate(gene_id_all = gene_id)) %>% 
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # take genes that were not included in the lastz filtered alignment 
    filter(lastzCov == "NA") %>%
    # tail()
    distinct(gene_id) %>%
    # nrow()
    write.table(paste0(venndiagOut, "EV_mazia_panMA.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  
  # Against MB reads mapping and/or CDS alignment
  read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"mazia_specific_against_panmb_nucmer_lastz.MB.CDS.bed"),
             header = F) %>%
    rename(
      seqid = V1,
      start = V2,
      end = V3,
      lastzCov = V4
    ) %>%
    distinct() %>%
    # nrow()
    # # join with CDS gene IDs and gene ontology (GO) terms. 
    # full_join(
    #   EV_MAB_reads_coverage_geneID_GO_terms %>% 
    #     filter(ref_genome == "ev_mazia") %>% 
    #     select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
    #     distinct() 
    # ) %>%
    # mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # filter(lastzCov == "NA") %>%
    # distinct(gene_id) %>%
    left_join(
      # Join with list of mazia specific genes from panma reads
      read.delim(paste0(result_dir,"filtered_ev_mazia_cds_less_25_per_panmb_cov.bed"), header = F) %>%
        mutate(source = "panMA_reads") %>%
        rename(seqid = V1,
               start = V2,
               end = V3)) %>%
    mutate(source = str_replace_na(source, "NA")) %>%
    # join with CDS gene IDs and gene ontology (GO) terms. 
    full_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_mazia") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov < 24.4) %>%
    distinct(gene_id,lastzCov) %>%
    # distinct(gene_id) %>%
    # nrow()
    # Join all gene ids of ref genome 
    full_join(EV_MAB_reads_coverage_geneID_GO_terms %>% 
                filter(ref_genome == "ev_mazia") %>% 
                distinct(gene_id, ref_genome) %>%
                mutate(gene_id_all = gene_id)) %>% 
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    # take genes that were not included in the lastz filtered alignment 
    filter(lastzCov == "NA") %>%
    # tail()
    distinct(gene_id) %>%
    # nrow()
    write.table(paste0(venndiagOut, "EV_mazia_panMB.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  # Against EV reads mapping 
  EV_MAB_reads_coverage %>% 
    filter(reads_used != "qfiltered" ) %>%
    filter(reads_genome != "S-sub genome banana genotypes",
           reads_genome != "T-sub genome banana genotypes") %>%
    # distinct(query) 
    
    mutate(fraction_overlap = round(fraction_overlap,2)) %>%
    filter(ref_genome == "ev_mazia") %>%
    # distinct(query)
    filter(query == "panEV_reads") %>%
    select(seq.name,start,end,fraction_overlap) %>%
    group_by(seq.name,start,end) %>%
    summarise(fraction_overlap = max(fraction_overlap)) %>%
    ungroup() %>% 
    # 11 EV (Mazia) genes showed < 25% panEV reads coverage but +25 CDS alignment coverage against Bedadeti CDS hence included as shared EV genes  
    # filter(fraction_overlap >= 0.25) %>%
    distinct(seq.name, start, end) %>%
    rename(
      seqid = seq.name
    ) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "ev_mazia") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    distinct(gene_id) %>% 
    write.table(paste0(venndiagOut, "EV_mazia_panEV.gene.txt"), col.names = F, sep = '\t', quote = F, row.names = F)
  
  
  
  ## Generate venn diagrams
  
  Gene_list_VennDiag_PanEV_MAB <- list()
  
  # OG_list_VennDiag_EV_MAB[["EV (Bedadeti)"]] <- read.table("../monocot_blastp/advanced_orthofinder/OrthogroupsVennDiag/Ensete_ventricosum_bedadeti.OG.txt", header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(venndiagOut, "EV_mazia_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(venndiagOut, "EV_mazia_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(venndiagOut, "EV_mazia_panMB.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
  
  # Plot Ven diagra 
  # library(VennDiagram)
  # library(ggVennDiagram)
  
  venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
  venn_data_process <- process_data(venn_data)
  
  EV_mazia_vs_MAB_vendiag.v2 <-
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
  
  
# Collect gene sets
  # EV genes that do not have EV reads support  
  venn_data_process$regionData$item[[4]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>% 
    # head() %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
      # EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      #   select(-GO_terms) %>%
      #   filter(reads_used != "qfiltered")  %>%
      #   filter(gene_repeat == "gene")
    ) %>%
    distinct(gene_id) %>%
    nrow ()
  
  
  ## Collect gene sets 
  
  EV_mazia_specific.v2 <-
    venn_data_process$regionData$item[[1]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  EV_mazia_MAB_core.v2 <-
    venn_data_process$regionData$item[[7]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  EV_mazia_MA_specific.v2 <-
    venn_data_process$regionData$item[[4]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  EV_mazia_MB_specific.v2 <-
    venn_data_process$regionData$item[[5]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  
  EV_mazia_noEVReadSupport_Butsyntenic2MAB <-
  venn_data_process$regionData$item[[6]] %>% 
    as.data.frame() %>% 
    rename(gene_id=".") %>% 
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms
    ) %>%
    select(GO_terms,gene_id)  %>%
    distinct() 
  

## For Musa acuminata against EV reads mapping and/or CDS alignment
  
  
  # Against EV  reads mapping and/or CDS alignment
  read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"MA_specific_against_panevReads_nucmer_lastz.mazia.CDS.bed"),
             header = F) %>%
    rbind(paste0(SPECIES_SPECIFIC_CDS_DIR,"MA_specific_against_panevReads_nucmer_lastz.bedadeti.CDS.bed"),
          header = F) %>%
    rename(
      seqid = V1,
      start = V2,
      end = V3,
      lastzCov = V4
    ) %>%
    distinct() %>%
    # join with CDS gene IDs and gene ontology (GO) terms. 
    mutate(start = as.numeric(start),
           end = as.numeric(end),
           lastzCov = as.numeric(lastzCov)) %>%
    # nrow()
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov != "NA") %>%
    mutate(lastzCov = as.numeric(lastzCov)) %>%
    filter(lastzCov < 24.4) %>%
  full_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "musa_acuminata") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>% 
    mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
    filter(lastzCov == "NA") %>%
    distinct(gene_id) %>% 
    write.table(paste0(venndiagOut, "MA_panEV.gene.txt"), 
                col.names = F, 
                sep = '\t', 
                quote = F, 
                row.names = F)
  
  
  # Against EV reads mapping and/or CDS alignment
  EV_MAB_reads_coverage %>% 
    filter(reads_used != "qfiltered" ) %>%
    filter(reads_genome != "S-sub genome banana genotypes",
           reads_genome != "T-sub genome banana genotypes") %>%
    # distinct(query) 
    
    mutate(fraction_overlap = round(fraction_overlap,2)) %>%
    filter(ref_genome == "musa_acuminata") %>% 
    # distinct(query)
    filter(query == "panMA_reads") %>% 
    select(seq.name,start,end,fraction_overlap) %>%
    group_by(seq.name,start,end) %>%
    summarise(fraction_overlap = max(fraction_overlap)) %>%
    ungroup() %>% 
    # 11 EV (Mazia) genes showed < 25% panEV reads coverage but +25 CDS alignment coverage against Bedadeti CDS hence included as shared EV genes  
    # filter(fraction_overlap >= 0.25) %>%
    distinct(seq.name, start, end) %>%
    rename(
      seqid = seq.name
    ) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms %>% 
        filter(ref_genome == "musa_acuminata") %>% 
        select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
        distinct() 
    ) %>%
    distinct(gene_id) %>% 
    write.table(paste0(venndiagOut, "MA_panMA.gene.txt"), 
                col.names = F, 
                sep = '\t', 
                quote = F, 
                row.names = F)
  

  ## Generate venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(venndiagOut, "MA_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nA-genome"]] <- read.table(paste0(venndiagOut, "MA_panMA.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

MA_vs_EV_vendiag.v2 <-
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

# MA genes that do not have MA reads support  
venn_data_process$regionData$item[[3]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  distinct(gene_id) %>%
  nrow()



## collect gene sets

MA_specific.v2 <-
  venn_data_process$regionData$item[[2]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 


MA_EV_core.v2 <-
  venn_data_process$regionData$item[[3]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_wide
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 


## For Musa balbisiana

# Against EV  reads mapping and/or CDS alignment
read.delim(paste0(SPECIES_SPECIFIC_CDS_DIR,"MB_specific_against_panevReads_nucmer_lastz.mazia.CDS.bed"),
           header = F) %>%
  rbind(paste0(SPECIES_SPECIFIC_CDS_DIR,"MB_specific_against_panevReads_nucmer_lastz.bedadeti.CDS.bed"),
        header = F) %>%
  rename(
    seqid = V1,
    start = V2,
    end = V3,
    lastzCov = V4
  ) %>%
  distinct() %>%
  # nrow()
  # join with CDS gene IDs and gene ontology (GO) terms. 
  mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
  mutate(lastzCov = as.numeric(lastzCov)) %>%
  filter(lastzCov < 24.4) %>%
  # distinct(seqid,start,end) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  full_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      filter(ref_genome == "musa_balbisiana") %>% 
      select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
      distinct() 
  ) %>%
  mutate(lastzCov = str_replace_na(lastzCov,"NA")) %>%
  # mutate(lastzCov = as.numeric(lastzCov)) %>%
  filter(lastzCov == "NA") %>%
  # head()
  # summary()
  # summarize(max (lastzCov))
  # filter(lastzCov == "NA") %>%
  distinct(gene_id) %>% 
  # nrow()
  write.table(paste0(venndiagOut, "MB_panEV.gene.txt"), 
              col.names = F, 
              sep = '\t', 
              quote = F, 
              row.names = F)


# Against EV reads mapping and/or CDS alignment
EV_MAB_reads_coverage %>% 
  filter(reads_used != "qfiltered" ) %>%
  filter(reads_genome != "S-sub genome banana genotypes",
         reads_genome != "T-sub genome banana genotypes") %>%
  # distinct(query) 
  
  mutate(fraction_overlap = round(fraction_overlap,2)) %>%
  filter(ref_genome == "musa_balbisiana") %>% 
  # distinct(query)
  filter(query == "panMB_reads") %>% 
  select(seq.name,start,end,fraction_overlap) %>%
  group_by(seq.name,start,end) %>%
  summarise(fraction_overlap = max(fraction_overlap)) %>%
  ungroup() %>% 
  # 11 EV (Mazia) genes showed < 25% panEV reads coverage but +25 CDS alignment coverage against Bedadeti CDS hence included as shared EV genes  
  # filter(fraction_overlap >= 0.25) %>%
  distinct(seq.name, start, end) %>%
  rename(
    seqid = seq.name
  ) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>% 
      filter(ref_genome == "musa_balbisiana") %>% 
      select(seqid, start, end, gene_id, GO_terms, ref_genome) %>% 
      distinct() 
  ) %>%
  distinct(gene_id) %>% 
  write.table(paste0(venndiagOut, "MB_panMB.gene.txt"), 
              col.names = F, 
              sep = '\t', 
              quote = F, 
              row.names = F)

# Figure of main venn diagrams

Gene_list_VennDiag_PanEV_MAB <- list()

Gene_list_VennDiag_PanEV_MAB[["EV"]] <- read.table(paste0(venndiagOut, "MB_panEV.gene.txt"), header = F, stringsAsFactors = FALSE)$V1
Gene_list_VennDiag_PanEV_MAB[["Musa spp\nB-genome"]] <- read.table(paste0(venndiagOut, "MB_panMB.gene.txt"), header = F, stringsAsFactors = FALSE)$V1

venn_data <- Venn(Gene_list_VennDiag_PanEV_MAB)
venn_data_process <- process_data(venn_data)

MB_vs_EV_vendiag.v2 <-
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
venn_data_process$regionData$item[[2]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>%
  # head() %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  distinct(gene_id) %>%
  nrow()

MB_specific.v2 <-
  venn_data_process$regionData$item[[2]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  select(GO_terms,gene_id)  %>%
  distinct() 

MB_EV_core.v2 <-
  venn_data_process$regionData$item[[3]] %>% 
  as.data.frame() %>% 
  rename(gene_id=".") %>% 
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms
  ) %>%
  # distinct(gene_id) %>%
  # nrow()
  select(GO_terms,gene_id)  %>%
  distinct() 

library(patchwork)

EV_mazia_vs_MAB_vendiag.v2+EV_bedadeti_vs_MAB_vendiag.v2
ggsave(paste0(venndiagOut, "EV_mazia_bedadeti_vendiag.allMapped.synteny.02.png"), width=10, height=4)
ggsave(paste0(venndiagOut, "EV_mazia_bedadeti_vendiag.allMapped.synteny.02.pdf"), width=10, height=4)

MA_vs_EV_vendiag.v2+MB_vs_EV_vendiag.v2
ggsave(paste0(venndiagOut, "MA_MB_vendiag.allMapped.synteny.02.png"), width=6, height=4)
ggsave(paste0(venndiagOut, "MA_MB_vendiag.allMapped.synteny.02.pdf"), width=6, height=4)

(EV_mazia_vs_MAB_vendiag.v2/MA_vs_EV_vendiag.v2)
# /(MA_vs_EV_vendiag|MB_vs_EV_vendiag)

ggsave(paste0(venndiagOut, "EV_mazia_top_vs_MA_bottom_genes_vendiag.allMapped.synteny.02.pdf"), width=9, height=9)

(EV_bedadeti_vs_MAB_vendiag.v2/MB_vs_EV_vendiag.v2)
# /(MA_vs_EV_vendiag|MB_vs_EV_vendiag)
ggsave(paste0(venndiagOut, "EV_bedadeti_top_vs_MB_bottom_genes_vendiag.allMapped.synteny.02.pdf"), width=9, height=9)

# cluster EVMBA genes 

EV_MAB_custered_genes.v4 <-
  EV_bedadeti_specific.v2 %>%
  mutate(clust = "EV_bedadeti_specific") %>% 
  rbind(
    EV_bedadeti_MAB_core.v2 %>%
      mutate(clust = "EV_bedadeti_vs_MA_shared"),
    EV_bedadeti_MA_specific.v2 %>%
      mutate(clust = "EV_bedadeti_vs_MA_shared"),
    EV_bedadeti_MB_specific.v2 %>%
      mutate(clust = "EV_bedadeti_vs_MB_shared"),
    # EV_bedadeti_noEVReadSupport_Butsyntenic2MAB %>%
    #   mutate(clust = "EV_bedadeti_noEVReadSupport_Butsyntenic2MAB"),
    EV_mazia_specific.v2 %>%
      mutate(clust = "EV_mazia_specific"),
    EV_mazia_MAB_core.v2 %>%
      mutate(clust = "EV_mazia_vs_MAB_shared"),
    EV_mazia_MA_specific.v2 %>%
      mutate(clust = "EV_mazia_vs_MA_shared"),
    EV_mazia_MB_specific.v2 %>%
      mutate(clust = "EV_mazia_vs_MB_shared"),
    # EV_mazia_noEVReadSupport_Butsyntenic2MAB%>%
    #   mutate(clust = "EV_mazia_noEVReadSupport_Butsyntenic2MAB"),
    MA_specific.v2 %>%
      mutate(clust = "MA_specific"),
    MA_EV_core.v2 %>%
      mutate(clust = "MA_EV_shared"),
    # MA_noMAReadSupport_Butsyntenic2EV %>%
    #   mutate(clust = "MA_noMAReadSupport_Butsyntenic2EV"),
    MB_specific.v2 %>%
      mutate(clust = "MB_specific"),
    MB_EV_core.v2 %>%
      mutate(clust = "MB_EV_shared")
    # MB_noMBReadSupport_Butsyntenic2EV %>%
    #   mutate(clust = "MB_noMBReadSupport_Butsyntenic2EV")
  ) %>%
  distinct(gene_id,clust) 


EV_MAB_custered_genes.v4 %>%
  group_by(clust) %>%
  summarise(count = n())
# ggsave(paste0(DATA_DIR, "plot_out/EV_bedadeti_top_vs_MB_bottom_genes_vendiag.allMapped.synteny.02.pdf"), width=9, height=9)

# save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.RData"))
# load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.RData"))
# 
# save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.01.RData"))
# # load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.01.RData"))
# load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.02.RData"))
# 
# load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.03.RData"))
# save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.03.RData"))

# save.image(paste0(DATA_DIR,"panEV_MAB_reads_coverage.06.RData"))
# load(paste0(DATA_DIR,"panEV_MAB_reads_coverage.06.RData"))
