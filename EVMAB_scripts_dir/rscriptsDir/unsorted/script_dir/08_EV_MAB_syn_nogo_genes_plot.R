

install.packages("rcartocolor")
library(rcartocolor)

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits <-
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  # select(-isomer) %>%
  mutate(
    genome = str_replace(genome, "MA","*M. acuminata*"),
    genome = str_replace(genome, "MB","*M. balbisiana*"),
    genome = str_replace(genome,"EV_mazia", "*E. ventricosum* (Mazia)"),
    genome = str_replace(genome,"EV_bedadeti", "*E. ventricosum* (Bedadeti)") ) %>%
  rename(description=stitle_bp_annot)


EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%distinct(genome) %>%dput()

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits$genome = factor(EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits$genome,
                                                                      levels = c("*E. ventricosum* (Mazia)", 
                                                                                 "*E. ventricosum* (Bedadeti)", 
                                                                                 "*M. acuminata*", 
                                                                                 "*M. balbisiana*"))


# blast annotation enset specific genes that lack GO-terms annotation 
read.delim(paste0(DATA_DIR,"/EV_MAB.nogo_noreads.ncbi_nr.diamond.blastp.e05.out.gz")) %>%
  rbind(
  read.delim(paste0(DATA_DIR,"/EV_MAB.nogo_noreads.uniprot_sprot.diamond.blastp.e05.out.gz")),
  read.delim(paste0(DATA_DIR,"/EV_MAB.nogo_noreads.atha_osa.diamond.blastp.e05.out.gz")))%>% 
  filter(str_detect(qseqid,"EVBD.1.033364")|
           str_detect(qseqid,"EVBD.1.051740")|
           str_detect(qseqid,"EVMZ.1.034720")) %>%
  # filter(str_detect(stitle,"hypersensitive")) 
  select(qseqid,qlen,slen,pident,length,stitle) %>%
  write.table(paste0(DATA_DIR,"EV_specific_no_blastp_homology.txt"),
              col.names = T, 
              row.names = F,
              quote = F,
              sep = "\t")


  


  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>% 
    mutate(
      description = str_replace(description, "hypersensitive induced reaction","Hypersensitive-induced response protein"),
      description = str_replace(description, "disease resistance protein","Disease resistance protein")) %>%
    
    filter(str_detect(gene_id,"EVBD.1.033364")|
             str_detect(gene_id,"EVBD.1.051740")|
             str_detect(gene_id,"EVMZ.1.034720")) %>%
    distinct() %>%
    group_by(genome,description) %>%
    summarise(gene_count = n()) %>%
    as.data.frame() 
    
    
  
  

  
  
  
  
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot <- 
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
    mutate(
      description = str_replace(description, "hypersensitive induced reaction","Hypersensitive-induced response protein"),
      description = str_replace(description, "disease resistance protein","Disease resistance protein")) %>%

    filter(str_detect(gene_id,"EVBD.1.033364")|
             # str_detect(gene_id,"EVBD.1.051740")|
             str_detect(gene_id,"EVMZ.1.034720")) %>%
    group_by(genome,description) %>%
    summarise(gene_count = n()) %>%
    rbind(
  
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
    mutate(
      description = str_remove(description,"VND-interacting 1,\\s+"),
      description = str_replace(description,"expressed protein","Expressed protein"),
      description = str_replace(description,"polyproteins", "Transposable elements (TE)/TE domains"),
      description = str_replace(description,"Transposable elements \\(TE\\)\\/TE domains", "Transposable elements (TE) or<br> TE domains encoding proteins"),
      description = str_replace(description,"Probable serine/threonine-protein kinase kinX", "Serine/threonine-protein<br>kinase kinX"),
      description = case_when(str_detect(description,"CASEIN LYTIC PROTEINASE B-P") ~ "CASEIN LYTIC PROTEINASE B-P",
                              TRUE ~ description)) %>%

    
    arrange(desc(description)) %>%
    filter(!str_detect(description,"mucin"),
           # !str_detect(description,"mucin-2-like")|
           !str_detect(description,"extensin-like"),
           !str_detect(description,"Periaxin"),
           !str_detect(description,"ENODL1"),
           !str_detect(description,"Nucleoporin NSP1")) %>%
    group_by(genome,description) %>%
    summarise(gene_count = n()) %>%
    filter(gene_count > 1)) %>%
    mutate(
      description = case_when(str_detect(description,"induced") ~ "Hypersensitive-induced<br>response protein",
                              TRUE ~ description),
    ) %>%

    filter(!str_detect(description,"Expressed")) 
  
  
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot %>%
    as.data.frame()
    
  
  # Arrange orders of genomes for plotting 
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot$genome = factor(EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot$genome,
                                                             levels = c("*E. ventricosum* (Mazia)", 
                                                                        "*E. ventricosum* (Bedadeti)", 
                                                                        "*M. acuminata*"))
  
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot %>%
    filter(!str_detect(genome,"acuminata")) %>%
    distinct() %>%
    
  
  ggplot() +
  geom_col(
    aes( x=reorder(description, gene_count), y=gene_count, fill = genome),
    position = "dodge"
  ) + 
  scale_y_continuous(breaks =c(1, 3,6,9,13))+
  scale_fill_carto_d(type = "qualitative", palette = "Vivid")+
  labs(fill="Genomic reads source",
       y="Number of genes") +
  coord_flip()+
  theme(
    axis.title.x = element_markdown(size = 10, face = "bold"),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_markdown(size = 10, face = "bold"),
    axis.text.x = element_markdown(size = 10,face = "bold"),
    # axis.title.x
    # panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_blank(),
    legend.key.size = unit(0.5,"cm"),
    # legend.spacing = unit(2,"cm"),
    # legend.direction = "vertical",
    legend.key.height = unit(0.6,"cm"),
    legend.box.spacing = unit(0,"cm"),
    # legend.box.margin = margin(-10,-10,-10,-10),
    legend.margin = margin(0,0,0,0),
    legend.text = element_markdown(size =9,face = "bold"),
    legend.position = "top"
    
  )
  
  ggsave(paste0(DATA_DIR, "plot_out/EV_nogo_genes_annotation.01.pdf"), width=5, height=3)
  ggsave(paste0(DATA_DIR, "plot_out/EV_nogo_genes_annotation.png"), width=6, height=4)
  
  
  # Musa acuminata
  
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits_plot %>%
    filter(str_detect(genome,"acuminata")) %>%
    ggplot() +
    geom_col(
      aes( x=reorder(description, gene_count), y=gene_count, fill = genome),
      position = "dodge"
    ) + 
    scale_y_continuous(breaks =c(1,2,3,4))+
    scale_fill_carto_d(type = "qualitative", palette = "Earth")+
    labs(fill="Genomic reads source",
         y="Number of genes") +
    coord_flip()+
    theme(
      axis.title.x = element_markdown(size = 10, face = "bold"),
      axis.title.y = element_blank(),
      axis.line.x = element_line(),
      axis.line.y = element_line(),
      axis.text.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(size = 10,face = "bold"),
      panel.background = element_rect(fill="#FFFFFF", color=NA),
      legend.title = element_blank(),
      legend.key.size = unit(0.5,"cm"),
      # legend.spacing = unit(2,"cm"),
      # legend.direction = "vertical",
      legend.key.height = unit(0.6,"cm"),
      legend.box.spacing = unit(0,"cm"),
      # legend.box.margin = margin(-10,-10,-10,-10),
      legend.margin = margin(0,0,0,0),
      legend.text = element_markdown(size =9,face = "bold"),
      legend.position = "top"
      
    )
  
  
  ggsave(paste0(DATA_DIR, "plot_out/MA_nogo_genes_annotation.02.pdf"), width=4, height=2)
  ggsave(paste0(DATA_DIR, "plot_out/MA_nogo_genes_annotation.png"), width=4, height=3)
  

  