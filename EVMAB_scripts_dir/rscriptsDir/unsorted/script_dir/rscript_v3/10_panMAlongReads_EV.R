# Assessing EV specific genes using long read genome sequences from Musa species 
# Data from 00_EV_MAB

# panMALongReads_EVMZgenome_cov
# panMALongReads_EVBDgenome_cov


panMALongReads_EVMZgenome_cov %>%
  rbind(panMBLongReads_EVMZgenome_cov) %>%
  select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
  rename(frac_ovlp= fraction_overlap,
         lndrc = landrace) %>%
  head()

    
panMALongReads_EVBDgenome_cov %>%
  rbind(panMBLongReads_EVBDgenome_cov) %>%
      select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
      rename(frac_ovlp= fraction_overlap,
             lndrc = landrace) %>%
  head()

# checking genes in MA/MB who lacks sufficient NGS reads support using self mapping of long reads mapping
panMALongReads_MAgenome_cov %>%
  head()

panMALongReads_MAgenome_cov %>%
  rename(seqid = seq.name) %>%
  select(seqid,
         start,
         end,
         fraction_overlap,
         ref_genome,
         landrace,
         reads_used) %>%
  # filter(str_detect(end,"17649118")) %>%
  left_join(

    EV_MAB_custered_genes %>%
      filter(clust == "MA_noMAReadSupport_Butsyntenic2EV") %>%
      left_join(
        EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
          filter(gene_repeat == "gene") %>%
          select(seqid, start, end, gene_id, fraction_overlap,gene_repeat) %>%
          distinct()
      ) %>%
      select(seqid, start, end, gene_id, gene_repeat,fraction_overlap) %>%
      filter(fraction_overlap < 0.25) %>%
      rename(frac_ovrNGS = fraction_overlap))%>%
  filter(gene_repeat == "gene") %>%
  filter(fraction_overlap < 0.25)
  tail()

  # Assessing EV specific genes who do not have sufficient reads support from pan-MAB NGS using long reawds of Musa spp. genotypes
  
  # Bedadeti 
  EV_MAB_custered_genes %>%
    filter(clust == "EV_bedadeti_specific") %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
        filter(gene_repeat == "gene") %>%
        select(seqid, start, end, gene_id, fraction_overlap,gene_repeat) %>%
        distinct()) %>%
    select(seqid, start, end, gene_id, gene_repeat,fraction_overlap) %>%
    filter(fraction_overlap < 0.25) %>%
    rename(frac_ovrNGS = fraction_overlap)%>%
  # head() %>%
    left_join(
      panMALongReads_EVBDgenome_cov %>%
        rbind(panMBLongReads_EVBDgenome_cov) %>%
        select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
        rename(seqid= seq.name) 
    ) %>%
    filter(fraction_overlap >= 0.25 ) %>%
    distinct(gene_id) %>%
    nrow()
  
        
  # Mazia 
  
  EV_MAB_custered_genes %>%
    filter(clust == "EV_mazia_specific") %>%
    # distinct(gene_id) %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
        filter(gene_repeat == "gene") %>%
        select(seqid, start, end, gene_id, fraction_overlap,gene_repeat) %>%
        distinct()) %>%
    select(seqid, start, end, gene_id, gene_repeat,fraction_overlap) %>%
    filter(fraction_overlap < 0.25) %>%
    rename(frac_ovrNGS = fraction_overlap)%>%
    # head() %>%
    left_join(
      panMALongReads_EVMZgenome_cov %>%
        rbind(panMBLongReads_EVMZgenome_cov) %>%
        select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
        rename(seqid= seq.name)) %>%
    filter(fraction_overlap >= 0.25 ) %>%
    # head()
    distinct(gene_id) %>%
    nrow()
  
  
  
  panMALongReads_MAgenome_cov %>%
    rename(seqid = seq.name) %>%
    select(seqid,
           start,
           end,
           fraction_overlap,
           ref_genome,
           landrace,
           reads_used) %>%
    # filter(str_detect(end,"17649118")) %>%
    left_join(
      
      EV_MAB_custered_genes %>%
        filter(clust == "MA_noMAReadSupport_Butsyntenic2EV") %>%
        left_join(
          EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
            filter(gene_repeat == "gene") %>%
            select(seqid, start, end, gene_id, fraction_overlap,gene_repeat) %>%
            distinct()
        ) %>%
        select(seqid, start, end, gene_id, gene_repeat,fraction_overlap) %>%
        filter(fraction_overlap < 0.25) %>%
        rename(frac_ovrNGS = fraction_overlap))%>%
    filter(gene_repeat == "gene") %>%
    filter(fraction_overlap < 0.25) 

# coverage using using panMAB reads

## EVMB specific genes using panEV and panMAB NGS reads 
EV_MAB_custered_genes %>%
  filter(clust == "EV_bedadeti_specific" |
           clust == "EV_mazia_specific" |
           clust == "MA_specific" |
           clust == "MB_specific") %>%
  head() %>%
  
  # Join NGS coverage details 
    left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_synteny %>%
      filter(mapping_reads_source == "cross_species_reads") 
  ) %>%
  
  # Join coverage estimates of EV genes using longreads of MS, MT and MA genomes
  
  left_join(
    panMALongReads_EVMZgenome_cov %>%
      select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
      rename(frac_ovlp= fraction_overlap,
             lndrc = landrace) %>%
      rbind(
        panMALongReads_EVBDgenome_cov %>%
          select(seq.name,start,end,fraction_overlap, ref_genome, landrace) %>%
          rename(frac_ovlp= fraction_overlap,
                 lndrc = landrace)) %>%
      filter(str_detect(lndrc,"MA_"))
    
  ) %>%
  distinct()




EV_MAB_reads_coverage_geneID_GO_terms_wide %>%
  head()
