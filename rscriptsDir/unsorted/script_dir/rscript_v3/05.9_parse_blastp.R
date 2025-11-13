parse_blast_outformat6 <- function(path_dir="files_directory",
                                   # file_prefix="files_prefix",
                                   file_ending_pattern="files ending pattern including the last separator"#,
                                   # strings_to_remove="strings to remove after the prefix"
                                   ){
  # lists of blastn files
  
  # path_dir <- result_dir
  # Names = list.files(path = result_dir, pattern = "uniprot.blastp.e05.out.gz") 

  Names = list.files(path = path_dir, pattern=file_ending_pattern) #%>% 
    # str_remove(strings_to_remove)
  
  # blastp
  
  blast_out = c()
  
  
  # i="EV_bedadeti_go_terms_biological_process.uniprot.blastp.e05.out.gz"
  for (i in Names) {
    

    BASENAME <- print(i) %>%
       str_remove("_go_terms_biological_process.uniprot.blastp.e05.out.gz")
    
    blast_out <-
      read.delim(paste0(path_dir,"/", i), header = F) %>%
      dplyr::rename (
        qseqid = V1,
        sseqid = V2,
        qlen = V3,
        slen = V4,
        qstart = V5,
        qend = V6,
        sstart = V7,
        send = V8,
        stitle = V9,
        pident = V10,
        length = V11,
        evalue = V12,
        bitscore = V13
      ) %>%
      mutate(source = BASENAME ,
             query_cov = round(length/qlen,1),
             sub_cov = round(length/slen,1)) %>%
      rbind(blast_out)
    
    # separate (quer_sub, into = c("query_genome", "subject_genome"), sep="\\." )
    
    # blast_out[[i]] <- blastp # add it to your list
    
  }
  
  
  # Unlist  
# 
#   blast_out =
#     do.call(rbind, blast_out)
#   
  
  # remove row.names()
  
  row.names(blast_out) <- NULL
  #rm(blastp)
  
  # save output as r object 
  
  return(blast_out)
}


list.files(path = result_dir, pattern = "uniprot.blastp.e05.out.gz") 

blast_out <- parse_blast_outformat6(path_dir = result_dir,
                                    file_ending_pattern = ".uniprot.blastp.e05.out.gz"#,
                                    # strings_to_remove = ".uniprot.blastp.e05.out.gz"
                                    )

blast_out.1 <- 
blast_out %>% 
  distinct() %>%
  # filter(str_detect(source,"EV_mazia")) %>%
  # head()
  mutate(
    stitle_bp_annot = stitle,
    stitle_bp_annot = str_replace(stitle_bp_annot, "isoform \\w+\\d+", "isoform"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+\\|\\w+\\|\\w+\\s+"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"\\s+\\w+=\\w+\\s+.*"),
    stitle_bp_annot = case_when(
      str_detect(stitle_bp_annot, "\\w+-rich protein") ~ "proline-rich protein",
      str_detect(stitle_bp_annot, "\\w+agnesium-dependent") ~ "magnesium-dependent phosphatase 1",
      str_detect(stitle_bp_annot, "26\\w prot") ~ "26S proteasome regulatory subunit",
      str_detect(stitle_bp_annot, "resistance") ~ "disease resistance protein",
      str_detect(stitle_bp_annot, "\\w+-H2 finger protein") ~ "RING-H2 finger protein",
      str_detect(stitle_bp_annot, "\\w+-\\w+ protein kinase kinase kinase") |
        str_detect(stitle_bp_annot, "\\w+-\\w+ kinase kinase kinase")  ~ "MAPKKK",
      str_detect(stitle_bp_annot, "\\w+-\\w+ protein kinase") ~ "MAPK",
      str_detect(stitle, "MAPKKK17") ~ "MAPKKK",
      
      
      str_detect(stitle_bp_annot, "\\w+-like arabinogalactan") ~ "fasciclin-like arabinogalactan proteins",
      str_detect(stitle_bp_annot, "\\winc finger") ~ "Zinc finger proteins",
      str_detect(stitle_bp_annot, "repetitive proline-rich cell wall protein")  ~ "repetitive proline-rich cell wall protein",
      str_detect(stitle_bp_annot, "\\w+ dehalogenase-like") ~ "Haloacid dehalogenase-like hydrolase",
      str_detect(stitle_bp_annot, "embryogenesis") ~ "Late embryogenesis abundant proteins",
      str_detect(stitle_bp_annot, "ARM") |
        str_detect(stitle_bp_annot, "arm") |
        str_detect(stitle_bp_annot, "Arm") ~ "ARM repeat protein interacting WITH ABF2 protein",
      str_detect(stitle_bp_annot, "\\w+pimelate epimerase") ~ "diaminopimelate epimerase",
      str_detect(stitle_bp_annot, "\\w+rredoxin") ~ "ferredoxin",
      str_detect(stitle_bp_annot, "peroxiredoxin") ~ "peroxiredoxins",
      str_detect(stitle, "\\w+tein kinase-like") ~ "protein kinase-like domains",
      str_detect(stitle, "Polyprotein") |
        str_detect(stitle, "polyprotein") ~ "polyproteins",
      str_detect(stitle_bp_annot, "\\w3 ubiquitin") ~ "E3 ubiquitin-protein ligase",
      str_detect(stitle_bp_annot,"\\w+ypothetical.*") ~ "Hypothetical proteins",
      str_detect(stitle_bp_annot, "\\w+ncharacterized protein.*") ~ "Uncharacterized proteins",
      str_detect(stitle_bp_annot, "\\w+named protein product") ~ "Unnamed protein product",
      str_detect( stitle_bp_annot,"Repetitive proline-rich cell wall protein") ~ "Repetitive proline-rich cell wall protein",
      
      str_detect(stitle_bp_annot, "\\w+asparian strip.*")  ~ "Casparian strip membrane protein",
      str_detect(stitle_bp_annot, "NDR1/HIN1-like protein")  ~ "NDR1/HIN1-like protein",
      str_detect(stitle_bp_annot, "1-aminocyclopropane-1-carboxylate oxidase homolog") ~ "1-aminocyclopropane-1-carboxylate oxidase homolog",
      str_detect(stitle_bp_annot, "Proteasome subunit beta") ~ "Proteasome subunit beta",
      str_detect(stitle_bp_annot, "Pathogenesis-related thaumatin-like protein") ~ "Pathogenesis-related thaumatin-like protein",
      
      str_detect(stitle_bp_annot, "Glucan endo-1,3-beta-glucosidase") ~ "Glucan endo-1,3-beta-glucosidase",
      str_detect(stitle_bp_annot, "biotin carboxyl carrier protein of acetyl-CoA.*") ~ "biotin carboxyl carrier protein of acetyl-CoA carboxylase 1, ",
      str_detect(stitle_bp_annot, "Thioredoxin") ~ "Thioredoxin/Thioredoxin-like proteins",
      str_detect(stitle_bp_annot, "AT-rich interactive domain-containing protein") ~ "AT-rich interactive domain-containing protein",
      str_detect(stitle_bp_annot, "RING-H2 finger protein") ~ "RING-H2 finger protein",
      str_detect(stitle_bp_annot, "Putative cell agglutination protein") ~ "Putative cell agglutination protein",
      str_detect(stitle_bp_annot, "\\w+aumatin")  |
        str_detect(stitle_bp_annot, "THAUMATIN") ~ "Thaumatin family profile",
      str_detect(stitle_bp_annot, "\\w+isease resistance") ~ "Disease resistance proteins (RPP8-like/RGAs)",
      str_detect(stitle_bp_annot, "Stigma-specific STIG1-like protein") ~ "Stigma-specific STIG1-like protein",
      str_detect(stitle_bp_annot, "Biotin carboxyl carrier protein of acetyl-CoA carboxylase") ~ "Biotin carboxyl carrier protein of acetyl-CoA carboxylase",
      str_detect(stitle_bp_annot, "heat shock protein") ~ "Heat shock protein domains",
      
      str_detect(stitle_bp_annot, "Myb-related") |
        str_detect(stitle_bp_annot, "MYB") |
        str_detect(stitle_bp_annot, "MYB-like transcription factor") |
        str_detect(stitle_bp_annot, "MYB31 transcription") ~ "MYB-related Transcription factors",
      
      str_detect(stitle_bp_annot, "Riboflavin_synthase_like") |
        str_detect(stitle_bp_annot, "ribE: riboflavin synthase, alpha subunit") |
        str_detect(stitle_bp_annot, "Riboflavin_synthase_like") |
        str_detect(stitle_bp_annot, "RIBOFLAVIN SYNTHASE") |
        # str_detect(stitle_bp_annot, "RIBOFLAVIN SYNTHASE") |
        str_detect(stitle_bp_annot, "Riboflavin synthase") ~ "Riboflavin synthase domain-like",
      
      str_detect(stitle_bp_annot, "Universal stress protein")  ~ "Universal stress proteins",
      str_detect(stitle_bp_annot, "ranscriptional regulator")  ~ "Transcriptional regulators (TAC1/EARS)",
      str_detect(stitle_bp_annot, "Alpha/beta hydrolase domain-containing protein")  ~ "Alpha/beta hydrolase domain-containing protein",
      str_detect(stitle_bp_annot, "Subtilisin-like protease")  ~ "Subtilisin-like proteases",
      str_detect(stitle_bp_annot, "Dolichyl-phosphate beta-glucosyltransferase")  ~ "Dolichyl-phosphate beta-glucosyltransferases",
      str_detect(stitle_bp_annot, "\\w+itochondrial import")  ~ "Mitochondrial import inner membrane translocase",
      str_detect(stitle_bp_annot, "\\w+inc-finger homeodomain.*")  ~ "Zinc-finger homeodomain protein",
      str_detect(stitle_bp_annot, "DNA-directed RNA polymera")  ~ "DNA-directed RNA polymerase II",
      str_detect(stitle_bp_annot, "RNA-directed DNA polymera")  ~ "RNA-directed DNA polymerase",
      str_detect(stitle_bp_annot, "keratin-associated protein")  ~ "keratin-associated proteins",
      str_detect(stitle_bp_annot, "mRNA turnover protein.*")  ~ "mRNA turnover proteins",
      str_detect(stitle_bp_annot, "zf-CCHC domain-containing protein.*")  ~ "zf-CCHC domain-containing proteins",
      str_detect(stitle_bp_annot, "wall-associated receptor kinase.*")  ~ "wall-associated receptor kinases",
      str_detect(stitle_bp_annot, "RPM1-interacting protein 4")  ~ "RPM1-interacting protein 4",
      str_detect(stitle_bp_annot, "\\w+ound-.*")  ~ "Wound-inducible basic proteins",
      str_detect(stitle_bp_annot, "U-box domain-containing")  ~ "U-box domain-containing protein",
      str_detect(stitle_bp_annot, "leucine-rich repeat")  ~ "leucine-rich repeat-containing protein / serine/threonine protein",
      str_detect(stitle_bp_annot, "OXIDOREDUCTASE, 2OG-FE II  OXYGENASE FAMILY PROTEIN") ~ "2OG-Fe(II) oxygenase superfamily",
      str_detect(stitle_bp_annot, "DNA \\(cytosine-5\\)-methyltransferase \\w+")  ~ "DNA (cytosine-5)-methyltransferase (DRMB/DRM1/DRM2/DRM3)",
      str_detect(stitle_bp_annot, "Sulfate/thiosulfate import ATP-binding protein")  ~ "Sulfate/thiosulfate import ATP-binding protein",
      str_detect(stitle_bp_annot, "PHD finger protein") ~ "PHD finger proteins (MALE STERILITY 1/MEIOCYTE DEATH 1)",
      
      str_detect(stitle_bp_annot, "\\(S\\)-8-oxocitronellyl enol synthase")  ~ "(S)-8-oxocitronellyl enol synthases",
      str_detect(stitle_bp_annot, "Zinc finger protein") ~  "Zink finger proteins",
      str_detect(stitle_bp_annot, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
      str_detect(stitle_bp_annot, "PLASMODESMATA CALLOSE-BINDING PROTEIN") ~  "PLASMODESMATA CALLOSE-BINDING PROTEIN",
      str_detect(stitle_bp_annot, "\\w+ranscription fact.*") ~  "Transcription factors",
      str_detect(stitle_bp_annot, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
      str_detect(stitle_bp_annot, "1-aminocyclopropane-1-carboxylate.*") ~  "1-aminocyclopropane-1-carboxylate",
      str_detect(stitle_bp_annot, "60S ribosomal.*") ~  "60S ribosomal proteins",
      # 
      str_detect(stitle_bp_annot, "Retrovirus-related") |
        str_detect(stitle_bp_annot, "\\w+ransposon") |
        str_detect(stitle_bp_annot, "Enzymatic polyprotein") |
        str_detect(stitle_bp_annot, "Polyprotein") |
        str_detect(stitle_bp_annot, "CHROMO DOMAIN-CONTAINING PROTEIN") |
        str_detect(stitle_bp_annot, "MuDR family transposase") |
        str_detect(stitle_bp_annot, "RNase_HI_RT_Ty3") |
        str_detect(stitle_bp_annot, "scriptase") |
        str_detect(stitle_bp_annot, "polymerases") |
        str_detect(stitle_bp_annot, "COPIA") |
        str_detect(stitle_bp_annot, "Genome polyprotein") |
        str_detect(stitle_bp_annot, "DNA BINDING PROTEIN") |
        str_detect(stitle_bp_annot, "POLYPROTEIN") |
        str_detect(stitle_bp_annot, "CHROMO") |
        str_detect(stitle_bp_annot, "Copia protein") |
        str_detect(stitle_bp_annot, "transposable") |
        str_detect(stitle_bp_annot, "Enzymatic polyprotein") |
        str_detect(stitle_bp_annot, "NAC domain-containing protein") |
        str_detect(stitle_bp_annot, "SPOSON") ~ "Transposable elements (TE)/TE domains",
      
      str_detect(stitle_bp_annot, "SMALL NUCLEAR RNA ACTIVATING") | 
        str_detect(stitle_bp_annot,"SNRNA-ACTIVATING PROTEIN") |
        str_detect(stitle_bp_annot,"SMALL NUCLEAR RNA ACTIVATING COMPLEX") |
        str_detect(stitle_bp_annot,"Small nuclear RNA activating complex") ~ "Small nuclear RNA activating complex",
      str_detect(stitle_bp_annot, "E3 UBIQUITIN")  |
        str_detect(stitle_bp_annot, "UBIQUITIN-PROTEIN LIGASE") |
        str_detect(stitle_bp_annot, "E3 UBIQUITIN-PROTEIN") ~  "E3 ubiquitin proteins ligases",
      str_detect(stitle_bp_annot, "Mitochondrial carrier") ~  "Mitochondrial carrier protein",
      str_detect(stitle_bp_annot, "\\w+acalin-type")| 
        str_detect(stitle_bp_annot, "\\w+acalin-related")~  "Jacalin-related lectin",
      str_detect(stitle_bp_annot, "Chaperone protein DnaJ") ~  "Chaperone protein DnaJ",
      
      str_detect(stitle_bp_annot, "\\w+rotein of unknown function") |
        str_detect(stitle_bp_annot, "Domain of unknown function") |
        str_detect(stitle_bp_annot, "UPF0481 protein") ~  "Protein of unknown function",
      str_detect(stitle_bp_annot, "5'-NUCLEOTIDASE DOMAIN-CONTAINING") ~  "5' nucleotidase family",
      # str_detect(stitle_bp_annot, "5'-NUCLEOTIDASE DOMAIN-CONTAINING") ~  "5' nucleotidase family",
      str_detect(stitle_bp_annot, "26S prot") ~  "26S proteasome regulatory",
      str_detect(stitle_bp_annot, "HMG boxes") |
        str_detect(stitle_bp_annot, "HMG-box domain") |
        str_detect(stitle_bp_annot, "HMG-box") ~ "High mobility group box domain",
      str_detect(stitle_bp_annot, "\\w+ucleoredoxin") ~  "Nucleoredoxin",
      str_detect(stitle_bp_annot, "NUCLEAR MIGRATION PROTEIN NUDC") ~  "NUCLEAR MOVEMENT PROTEIN NUDC",
      str_detect(stitle_bp_annot, "BSD domain") |
        str_detect(stitle_bp_annot, "BSD DOMAIN") ~  "BSD domain-like",
      str_detect(stitle_bp_annot, "ALUMINUM INDUCED PROTEIN") ~  "Aluminium induced protein",
      str_detect(stitle_bp_annot, "P-loop containing nucleoside triphosphate") ~  "P-loop containing nucleoside triphosphate hydrolase",
      str_detect(stitle_bp_annot, "TRANSLATION ELONGATION") ~  "Translational elongation factor related protein",
      str_detect(stitle_bp_annot, "ARM REPEAT SUPERFAMILY PROTEIN") |
        str_detect(stitle_bp_annot, "ARM repeat") ~  "ARM REPEAT SUPERFAMILY PROTEIN",
      str_detect(stitle_bp_annot, "CHROMO DOMAIN-CONTAINING PROTEIN") |
        str_detect(stitle_bp_annot, "Reverse transcriptase") |
        str_detect(stitle_bp_annot, "DNA-binding pseudobarrel domain") |
        str_detect(stitle_bp_annot, "Ribonuclease H") |
        str_detect(stitle_bp_annot, "retrotransposon") |
        str_detect(stitle_bp_annot, "reverse transcriptase") |
        str_detect(stitle_bp_annot, "retrotransposable") |
        str_detect(stitle_bp_annot, "Putative enzymatic polyprotein") |                      
        str_detect(stitle_bp_annot, "Integrase catalytic domain-containing protein") |
        str_detect(stitle_bp_annot, "gag/pol protein") |
        str_detect(stitle_bp_annot, "RNA-DIRECTED DNA POLYMERASE HOMOLOG") ~  "Transposable elements (TE)/TE domains" ,
      
      str_detect(stitle_bp_annot, "proline-rich") |
        str_detect(stitle_bp_annot, "basic proline-rich protein-like") |
        str_detect(stitle_bp_annot, " Proline rich extensin signature") ~  "Proline-rich protein",
      str_detect(stitle_bp_annot, "\\w+lutamine-rich protein") ~  "Glutamine-rich protein 2-like",
      
      str_detect(stitle_bp_annot, "Zink finger proteins") |
        str_detect(stitle_bp_annot, "Zinc finger CCHC-type superfamily") ~  "Zinc finger proteins",
      str_detect(stitle_bp_annot, "DNA-damage-repair.*") |
        str_detect(stitle_bp_annot, "DNA damage-repair.*") ~  "DNA-damage-repair/toleration protein",
      
      # str_detect(stitle_bp_annot, "ALUMINUM INDUCED PROTEIN") ~  "Aluminium induced protein",
      str_detect(stitle_bp_annot, "SCO\\d+-like protein") ~  "SCO1 or SCO2-like protein",
      
      str_detect(stitle_bp_annot, "DNA  CYTOSINE-5- -METHYLTRANSFERASE 3-RELATED") ~  "DNA (CYTOSINE-5)-METHYLTRANSFERASE DRM",
      str_detect(stitle_bp_annot, "Universal stress proteins") ~  "Universal stress protein family",
      # str_detect(stitle_bp_annot, "Universal stress proteins") ~  "Universal stress protein family",
      str_detect(stitle_bp_annot, "\\w+entapeptide") ~  "Pentapeptide repeat-containing proteins",
      str_detect(stitle_bp_annot, "\\w+arge ribosomal subunit") ~  "Large ribosomal subunit proteins",
      str_detect(stitle_bp_annot, "\\annose/glucose-specific") ~  "Mannose/glucose-specific lectin",
      str_detect(stitle_bp_annot, "coumaroyl CoA and feruloyl") ~  "Bi-functional coumaroyl CoA and feruloyl CoA ortho-hydroxylase",
      str_detect(stitle_bp_annot, "\\w+athogenesis-related protein")  |
        str_detect(stitle_bp_annot, "\\w+athogenesis-related 5 protein") ~  "Pathogenesis-related proteins",
      str_detect(stitle_bp_annot, "\\w+usceptibility protein") ~  "Diease susceptibility protein",
      str_detect(stitle_bp_annot, "\\w+esquiterpene synthase") ~  "Sesquiterpene synthase",
      str_detect(stitle_bp_annot, "SCO\\d+\\-like protein") ~  "SCO1 or SCO2-like protein",
      str_detect(stitle_bp_annot, "LRR receptor-like serine") ~  "LRR receptor-like serine/threonine-protein kinase",
      str_detect(stitle_bp_annot, "inorganic phosphate transporter") ~  "Inorganic phosphate transporter",
      str_detect(stitle_bp_annot, "Acyl-CoA-binding domain-containing protein") ~  "Acyl-CoA-binding domain-containing protein",
      str_detect(stitle_bp_annot, "Protein EXORDIUM-like") ~  "Protein EXORDIUM-like",
      str_detect(stitle_bp_annot, "CBL-interacting serine") ~  "CBL-interacting serine/threonine-protein kinase",
      str_detect(stitle_bp_annot, "\\w+oluble starch synthase") ~  "Soluble starch synthase",
      str_detect(stitle_bp_annot, "chain acyl-CoA synthetase") ~  "Long chain acyl-CoA synthetase",
      str_detect(stitle_bp_annot, "Long-chain-fatty-acid--CoA ligase") ~  "Long-chain-fatty-acid--CoA ligase",
      str_detect(stitle_bp_annot, "BTB/POZ domain and ankyrin") ~  "BTB/POZ domain and ankyrin repeat-containing protein",
      str_detect(stitle_bp_annot, "chloroplastic") ~  "Chloroplastic proteins",
      str_detect(stitle_bp_annot, "\\w+lutathione S-transferase") ~  "Glutathione S-transferase",
      str_detect(stitle_bp_annot, "Abscisic acid receptor") ~  "Abscisic acid receptor",
      str_detect(stitle_bp_annot, "Protein DETOXIFICATION") ~  "Protein DETOXIFICATION",
      str_detect(stitle_bp_annot, "globulin seed storage protein") ~  "11S or 13S globulin seed storage protein",
      str_detect(stitle_bp_annot, "ABC transporter G family") ~  "ABC transporter G family member",
      str_detect(stitle_bp_annot, "4-coumarate--CoA") ~  "4-coumarate--CoA ligase",
      str_detect(stitle_bp_annot, "\\w+ypersensitive-induced response protein") ~  "Hypersensitive-induced response protein",
      str_detect(stitle_bp_annot, "ATP-dependent RNA helicase") ~  "ATP-dependent RNA helicase",
      str_detect(stitle_bp_annot, "\\w+norganic phosphate transporter") ~  "Inorganic phosphate transporte",
      str_detect(stitle_bp_annot, "Jasmonate-induced oxygenas") ~  "Jasmonate-induced oxygenas",
      str_detect(stitle_bp_annot, "Protein FLOWERING")  |
        str_detect(stitle_bp_annot, "RICE FLOWERING") ~  "Protein FLOWERING LOCUS T",
      str_detect(stitle_bp_annot, "Protein VERNALIZATION") ~  "Protein VERNALIZATION",
      str_detect(stitle_bp_annot, "Fatty acyl-CoA synthetase") ~  "Fatty acyl-CoA synthetase",
      str_detect(stitle_bp_annot, "Feruloyl CoA ortho-hydroxyla") ~  "Feruloyl CoA ortho-hydroxyla",
      str_detect(stitle_bp_annot, "Protein DMR6-LIKE OXYGENASE") ~  "Protein DMR6-LIKE OXYGENASE",
      str_detect(stitle_bp_annot, "Short chain dehydrogenase") ~  "Short chain dehydrogenase",
      str_detect(stitle_bp_annot, "Glutelin type-\\w+") ~  "Glutelin type-A or B",
      str_detect(stitle_bp_annot, "Chaperone protein") ~  "Chaperone protein",
      str_detect(stitle_bp_annot, "RGF1 INDUCIBLE TRANSCRIPTION FACTOR") ~  "RGF1 INDUCIBLE TRANSCRIPTION FACTOR",
      str_detect(stitle_bp_annot, "\\s+\\d+") ~  "",
      # str_detect(stitle_bp_annot, "\\s+AFC") ~  "",
      # str_detect(stitle_bp_annot, "\\s+CLK\\s+") ~  "",
      # str_detect(stitle_bp_annot, "\\s+clk\ds+") ~  "",
      # str_detect(stitle_bp_annot, "\\s+\\d+") ~  "",
      str_detect(stitle_bp_annot, "\\s+\\w+\\d+") ~  "",
      str_detect(stitle_bp_annot, "\\s+\\w+.\\d+") ~  "",
      str_detect(stitle_bp_annot, "-\\w+\\d+") ~  "",
      # str_detect(stitle_bp_annot, "\\s+HAC") ~  "",
      # str_detect(stitle_bp_annot, "\\arge") ~  "Large",
      # str_detect(stitle_bp_annot, "\\arge") ~  "Large",
      # str_detect(stitle_bp_annot, "\\arge") ~  "Large",
      
      
      TRUE ~ stitle_bp_annot))%>%
  distinct(qseqid, stitle_bp_annot) %>%
  # group_by(qseqid, stitle_bp_annot) %>%
  as.data.frame() %>%
  filter(stitle_bp_annot != "") %>%
  separate(qseqid, into = c("seq.name", "isomer"), sep = "-") %>%
  # head()
  left_join(

interproscan_GO_terms %>%
  distinct() %>%
  filter(str_detect(GO_terms,"GO")) %>%
  # filter(lastz_cov < 0.25) %>%
  # summary
  left_join(
    read.delim(paste0(result_dir,"go_terms_table.txt"), header = T)) %>%
  filter(aspects == "biological_process") %>%
  filter(!str_detect(description,"obsolete")) %>%
  filter(clust == "MA_specific" |
           clust == "MB_specific" |
           clust == "EV_mazia_specific" |
           clust == "EV_bedadeti_specific") %>%
  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms %>%
      select(-GO_terms)
  ) %>%
  select(seqid, start, end, gene_id, clust, lastz_cov, GO_terms, description) %>%
  distinct() %>%
  # filter(#clust == "MA_specific" 
  #   #clust == "MB_specific" 
  #   clust == "EV_mazia_specific" 
  #     # clust == "EV_bedadeti_specific"
  # ) %>% 
  rename(seq.name = gene_id)) %>%
  # filter(#description=="carbohydrate metabolic process" |
  #          str_detect(description, "glycerol-3-phosphate") 
  #          # str_detect(description, "glycolytic process") 
  #          # str_detect(description, "defense")
  #          # str_detect(description, "glutathione metabolic process")
  #        )
  # distinct(stitle_bp_annot)
  # filter(
  #   #description=="carbohydrate metabolic process" |
  #   str_detect(stitle_bp_annot, "ransposable") |
  #   # str_detect(description, "ranscription")
  #   str_detect(description, "DNA integration")
  #   # str_detect(description, "glutathione metabolic process")
  # ) %>%
  # select(seqid, start, end,seq.name) %>%
  # filter(!str_detect(seq.name,"EVMZ")) %>%
  # write.table(paste0(result_dir,"MA_specific_DNA_intergration_transposable.cds.bed"),
  #             col.names = F,
  #             row.name = F,
  #             quote = F,
  #             sep = "\t")
  distinct()



  





