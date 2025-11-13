kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, bedadti_musa_ac_align_cov_grange_chr1,  
       y1=bedadti_musa_ac_align_cov_grange_chr1$id, 
       ymax=max(bedadti_musa_ac_align_cov_grange_chr1$id)/2, 
       col="#006400", r0=0.5, r1=1, border=NA)


kpBars(kp, centro_density,  
       y1=centro_density$id, 
       ymax=max(centro_density$id)/0.15, 
       col="red", r0=0.4, r1=1,border=NA)


kpBars(kp, bedadti_musa_ac_align_cov_grange_chr1,  
       y1=bedadti_musa_ac_align_cov_grange_chr1$id, 
       ymax=max(bedadti_musa_ac_align_cov_grange_chr1$id)/3.5, 
       col="#0000FF", r0=0.3, r1=0,border=NA)


centro_density <- 
  
  read.delim("C:/Users/Sadik/OneDrive - University of Warwick/Sadik/Ensete_banana_temp_dir/ALL_centro.density.100kb", header= F, sep = ' ') %>%
  dplyr::rename(chr = V1, 
                start = V2, 
                end = V3, 
                density = V4) %>%  
  dplyr::filter(str_detect(chr,"^A")) %>%
  mutate(chr = str_replace(chr, "A","chr")) %>% 
  filter(chr == "chr01") %>%
  bed_to_granges () 