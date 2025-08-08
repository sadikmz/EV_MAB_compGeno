
EV_MAB_reads_coverage_no_lastz_nucmer_alignment %>% 
  as.data.frame() %>%
  filter(
    group !="*Musa spp.* A-genome NGS reads vs <br> MA genome",
    group !="*Musa spp.* B-genome NGS vs <br> MB genome",
    group !="*E. ventricosum* NGS reads vs <br> MA genome",
    group !="*E. ventricosum* NGS reads vs <br> MB genome",
    gene_repeat=="gene"
  ) %>% 
  distinct() %>%
  filter(reads_used =="allMapped") %>%
  select( seq.name, start,  end, ref_genome, reads_used, gene_repeat, query, group, fraction_overlap, syn_cov) %>%
  mutate(syn_cov = str_replace_na(syn_cov,"NA")) %>%
  # filter(syn_cov=="NA") %>%
  filter(fraction_overlap < 0.245,
         syn_cov >=0.245


mutate(
    
    frac_ovrlp = fraction_overlap ,
    frac_ovrlp_updated = case_when(fraction_overlap < 0.245 & syn_cov >=0.245 ~ syn_cov,
                                   TRUE ~ frac_ovrlp),
    # frac_ovrlp_updated = case_when(fraction_overlap >= 0.245 ~ frac_ovrlp,
    #                                TRUE ~ frac_ovrlp_updated),
    frac_ovrlp_updated = str_replace_na(frac_ovrlp_updated,"NA")
  ) %>%
  filter(frac_ovrlp_updated=="NA") %>%
  head()

  
filter(fraction_overlap < 0.25) %>%
  head()
  head()
  group_by( seq.name, start,  end, ref_genome, reads_used, gene_repeat, query, group) %>%
  summarise(fraction_overlap = mean(fraction_overlap)) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot(aes(fraction_overlap)) +
  geom_freqpoly(linewidth=1, color="blue")+
  facet_wrap(~group, ncol = 2, nrow = 2, scales = "free")+
  labs(x = "Fraction coverage of each non-syntenic *E. ventricosum* gene<br>against genomic reads of either *Musa spp* (top) or *Ensete ventricosum* genotypes (bottom)" ,
       y = "Count of genes",
       color = "Source of EV<br>genomic reads (landraces)")+
  theme_bw()+
  theme(
    axis.title.x = element_markdown(size = 8,face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_markdown(size = 8,face ="bold"),
    axis.text.y = element_text(size = 8,face = "bold"),
    strip.text.x = element_markdown(size=8, face = "bold"),
    legend.text = element_markdown(margin = margin(r=8),size = 8,face = "bold"),
    legend.title = element_markdown(size=8, face="bold"),
    legend.key.size = unit(0.35,"cm"),
    legend.position = "none"
    
  )
  
  
  EV_MAB_genes %>%
    distinct() %>%
    group_by(genome) %>%
    summarise(count = n())
  