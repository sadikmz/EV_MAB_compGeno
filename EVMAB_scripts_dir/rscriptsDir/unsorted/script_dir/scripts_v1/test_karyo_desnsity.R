load("C:/Users/Sadik/OneDrive - University of Warwick/Sadik/Ensete_banana_temp_dir/karyo_data.RData")

kp <- plotKaryotype(genome=MA_genome_bed, ideogram.plotter = NULL, plot.type=6 )
# kpAddCytobandsAsLine(kp)
kpBars(kp, bedadti_musa_ac_align_cov_grange_chr1,
       y1=bedadti_musa_ac_align_cov_grange_chr1$id,
       ymax=max(bedadti_musa_ac_align_cov_grange_chr1$id)/2,
       col="#006400", r0=0.55, r1=1, border=NA)

kpBars(kp, centro_density,
       y1=centro_density$id,
       ymax=max(centro_density$id)/0.2,
       col="red", r0=0.42, r1=1,border=NA)

kpBars(kp, bedadti_musa_ac_align_cov_grange_chr1,
       y1=bedadti_musa_ac_align_cov_grange_chr1$id,
       ymax=max(bedadti_musa_ac_align_cov_grange_chr1$id)/3.2,
       col="#0000FF", r0=0.35, r1=0,border=NA)
