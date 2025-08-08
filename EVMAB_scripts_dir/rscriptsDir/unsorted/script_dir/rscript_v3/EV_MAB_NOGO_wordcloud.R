text <- readLines(file.choose())

# Read the text file from internet
# filePath <- "http://www.sthda.com/sthda/RDoc/example-files/martin-luther-king-i-have-a-dream-speech.txt"
# text <- readLines(filePath)

# text < Musa_AB_specifi_interproscan_NOGO_blastp$stitle

# MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro
# EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v1 <- 

read.delim("Z:sadik/evmab_dir/EV_MAB.nogo_noreads.ncbi_nr.diamond.blastp.e05.out.gz") %>%
# EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro %>%
  # mutate(stitle_sp = str_replace_na(stitle_sp, "NA")) %>%
  distinct() %>%
mutate(
  stitle_bp_annot = stitle,
  stitle_bp_annot = str_replace(stitle_bp_annot, "isoform \\w+\\d+", "isoform"),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+-rich protein") ~ "proline-rich protein",
    TRUE ~ stitle_bp_annot
  ),

  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+agnesium-dependent") ~ "magnesium-dependent phosphatase 1",
    TRUE ~ stitle_bp_annot
  ),

  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "26\\w prot") ~ "26S proteasome regulatory subunit",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "resistance") ~ "disease resistance protein",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+-H2 finger protein") ~ "RING-H2 finger protein",
    TRUE ~ stitle_bp_annot
  ),
  # RING-H2 finger protein
  
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+-\\w+ protein kinase kinase kinase") |
      str_detect(stitle_bp_annot, "\\w+-\\w+ kinase kinase kinase")  ~ "MAPKKK",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+-\\w+ protein kinase") ~ "MAPK",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle, "MAPKKK17") ~ "MAPKKK",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+-like arabinogalactan") ~ "fasciclin-like arabinogalactan proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\winc finger") ~ "Zinc finger proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "repetitive proline-rich cell wall protein")  ~ "repetitive proline-rich cell wall protein",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+ dehalogenase-like") ~ "Haloacid dehalogenase-like hydrolase",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "embryogenesis") ~ "Late embryogenesis abundant proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "ARM") |
      str_detect(stitle_bp_annot, "arm") |
      str_detect(stitle_bp_annot, "Arm") ~ "ARM repeat protein interacting WITH ABF2 protein",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+pimelate epimerase") ~ "diaminopimelate epimerase",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w+rredoxin") ~ "ferredoxin",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "peroxiredoxin") ~ "peroxiredoxins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle, "\\w+tein kinase-like") ~ "protein kinase-like domains",
    TRUE ~ stitle
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle, "Polyprotein") |
      str_detect(stitle, "polyprotein") ~ "polyproteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "\\w3 ubiquitin") ~ "E3 ubiquitin-protein ligase",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Solanum tuberosum\\]"
    ) ~ "Solanum tuberosum hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Solanum commersonii\\]"
    ) ~ "Solanum commersonii hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Citrus sinensis\\]"
    ) ~ "Citrus sinensis hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Ensete ventricosum\\]"
    ) ~ "Ensete ventricosum hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Shorea leprosula\\]"
    ) ~ "Shorea leprosula hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Prunus dulcis\\]") ~ "Prunus dulcis hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Musa balbisiana\\]"
    ) ~ "Musa balbisiana hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Malus baccata\\]") ~ "Malus baccata hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Actinidia rufa\\]") ~ "Actinidia rufa hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Cocos nucifera\\]") ~ "Cocos nucifera hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Trifolium subterraneum\\]"
    ) ~ "Trifolium subterraneum hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Prunus dulcis\\]") ~ "Prunus dulcis hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Arachis hypogaea\\]"
    ) ~ "Arachis hypogaea hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Musa balbisiana\\]"
    ) ~ "Musa balbisiana hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "hypothetical protein.*\\[Musa acuminata\\]") ~ "Musa acuminata hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),

  stitle_bp_annot = str_remove(stitle_bp_annot, "\\[.*"),
  stitle_bp_annot = str_remove(stitle_bp_annot, "PREDICTED:\\s+"),

  stitle_bp_annot = case_when(
    str_detect(
      stitle_bp_annot,
      "hypothetical protein.*\\[Mucuna pruriens\\]"
    ) ~ "Cocos nucifera hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "^hypothetical protein.*") ~ "hypothetical proteins",
    TRUE ~ stitle_bp_annot
  ),
  stitle_bp_annot = case_when(
    str_detect(stitle_bp_annot, "^uncharacterized protein.*") ~ "uncharacterized proteins",
    TRUE ~ stitle_bp_annot
  ),
  
  stitle_bp_annot = str_remove(stitle_bp_annot,"^\\w+.\\d+\\s+"),
  stitle_bp_annot = str_remove(stitle_bp_annot,"isoform.*"),
  stitle_bp_annot = str_remove(stitle_bp_annot,"LOC\\d.*+"),
  stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+_\\d+.*"),
  stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+_\\w+.*"),
  

  ) %>%
  distinct(stitle_bp_annot)
  
  
  ## interproscan_annotaiton 
  annotate_intpr = signature_description,
  annotate_intpr = case_when(
    str_detect(annotate_intpr, "^-") ~ interpro_description,
    str_detect(signature_accession, "PTHR24559:SF319") ~ "TRANSPOSON TY3-I GAG-POL POLYPROTEIN",
    str_detect(signature_accession, "PTHR46148") ~ "CHROMO DOMAIN-CONTAINING PROTEIN",
    # str_detect(annotate,"SK") ~ "Shikimate kinase",
    str_detect(signature_accession, "PTHR34629") ~ "PROLINE-RICH EXTENSIN-LIKE PROTEIN EPR1",
    str_detect(
      annotate_intpr,
      "P-loop containing nucleoside triphosphate hydrolase"
    ) ~ "P-loop containing nucleoside triphosphate hydrolases",

    str_detect(annotate_intpr, "\\wink \\winger") |
      str_detect(annotate_intpr, "ZINK FINGER") |
      str_detect(annotate_intpr, "ZINC FINGER PROTEIN CONSTANS-LIKE") |
      str_detect(annotate_intpr, "ZINC FINGER PROTEIN ZAT5") |
      str_detect(annotate_intpr, "FCS-LIKE ZINC FINGER 1-RELATED") |
      str_detect(annotate_intpr, "\\w+inc") ~ "Zink finger domains",
    
    str_detect(annotate_intpr, "heat shock protein") |
      str_detect(annotate_intpr, "17.8 kDa class I heat shock protein") ~ "Heat shock protein domains",
    str_detect(annotate_intpr, "\\waumatin")  |
      str_detect(annotate_intpr, "THAUMATIN") ~ "Thaumatin family profile",

    str_detect(annotate_intpr, "ASIC-LEUCINE ZIPPER") |
      str_detect(annotate_intpr, "\\yb") |
      str_detect(annotate_intpr, "MYB") ~ "Transcription factor families (BZIP/MYB)",
    
    str_detect(annotate_intpr, "Riboflavin_synthase_like") |
      str_detect(annotate_intpr, "ribE: riboflavin synthase, alpha subunit") |
      str_detect(annotate_intpr, "Riboflavin_synthase_like") |
      str_detect(annotate_intpr, "RIBOFLAVIN SYNTHASE") |
      # str_detect(annotate_intpr, "RIBOFLAVIN SYNTHASE") |
      str_detect(annotate_intpr, "Riboflavin synthase") ~ "Riboflavin synthase domain-like",
    
    str_detect(annotate_intpr, "UNIVERSAL STRESS PROTEIN PHOS34-LIKE") |
      str_detect(annotate_intpr, "USP_Like")  ~ "Universal stress protein family",
    
    str_detect(annotate_intpr, "OXIDOREDUCTASE, 2OG-FE II  OXYGENASE FAMILY PROTEIN") ~ "2OG-Fe(II) oxygenase superfamily",
    
    str_detect(annotate_intpr, "Thioredoxin") |
      str_detect(annotate_intpr, "THIOREDOXIN") |
      str_detect(annotate_intpr, "TRX_family") |
      str_detect(annotate_intpr, "TRX_family") |
      str_detect(annotate_intpr, "Thioredoxin-like") |
      str_detect(annotate_intpr, "THIOREDOXIN-LIKE PROTEIN") |
      str_detect(annotate_intpr, "USP_Like")  ~ "Universal stress protein family",
    
    str_detect(annotate_intpr, "sposon") |
      str_detect(annotate_intpr, "Reverse transcriptase/Diguanylate cyclase domain") |
      str_detect(annotate_intpr, "DNA-binding pseudobarrel domain superfamily") |
      str_detect(annotate_intpr, "DNA-binding pseudobarrel domain") |
      str_detect(annotate_intpr, "CHROMO DOMAIN-CONTAINING PROTEIN") |
      str_detect(annotate_intpr, "MuDR family transposase") |
      str_detect(annotate_intpr, "RNase_HI_RT_Ty3") |
      str_detect(annotate_intpr, "scriptase") |
      str_detect(annotate_intpr, "polymerases") |
      str_detect(annotate_intpr, "COPIA") |
      str_detect(annotate_intpr, "nuclease") |
      str_detect(annotate_intpr, "DNA BINDING PROTEIN") |
      str_detect(annotate_intpr, "POLYPROTEIN") |
      str_detect(annotate_intpr, "CHROMO") |
      str_detect(annotate_intpr, "copia") |
      str_detect(annotate_intpr, "TRANSPOSON") ~ "Transposable elements (TE)/TE domains",
    TRUE ~ annotate_intpr),
  annotate_intpr = str_remove(annotate_intpr, ";"),
  
  ## swiss-prot blastp against curated proteins 
  
  # sp|Q9FY82|NAC82_ARATH NAC domain-containing protein 82 OS=Arabidopsis thaliana OX=3702 GN=NAC082 PE=1 SV=1
  
  annotate_sp = stitle_sp,
  
  annotate_sp = str_remove(annotate_sp,"\\w+\\|\\w+\\|\\w+\\s+"),
  annotate_sp = str_remove(annotate_sp,"\\w+=\\w+\\s+.*"),
  annotate_sp = case_when(

    str_detect(
      annotate_sp,
      "Repetitive proline-rich cell wall protein"
    ) ~ "Repetitive proline-rich cell wall protein",

    str_detect(annotate_sp, "Casparian strip membrane protein")  ~ "Casparian strip membrane protein",
    str_detect(annotate_sp, "NDR1/HIN1-like protein")  ~ "NDR1/HIN1-like protein",
    str_detect(annotate_sp, "1-aminocyclopropane-1-carboxylate oxidase homolog") ~ "1-aminocyclopropane-1-carboxylate oxidase homolog",
    str_detect(annotate_sp, "Proteasome subunit beta") ~ "Proteasome subunit beta",
    str_detect(annotate_sp, "Pathogenesis-related thaumatin-like protein") ~ "Pathogenesis-related thaumatin-like protein",
    
    str_detect(annotate_sp, "Glucan endo-1,3-beta-glucosidase") ~ "Glucan endo-1,3-beta-glucosidase",
    # str_detect(annotate_sp, "heat shock protein") ~ "heat shock protein",
    str_detect(annotate_sp, "Thioredoxin") ~ "Thioredoxin/Thioredoxin-like proteins",
    str_detect(annotate_sp, "AT-rich interactive domain-containing protein") ~ "AT-rich interactive domain-containing protein",
    str_detect(annotate_sp, "RING-H2 finger protein") ~ "RING-H2 finger protein",
    str_detect(annotate_sp, "Putative cell agglutination protein") ~ "Putative cell agglutination protein",
    str_detect(annotate_sp, "\\waumatin")  |
      str_detect(annotate_sp, "THAUMATIN") ~ "Thaumatin family profile",
    
    
    str_detect(annotate_sp, "\\wisease resistance") ~ "Disease resistance proteins (RPP8-like/RGAs)",
    str_detect(annotate_sp, "Stigma-specific STIG1-like protein") ~ "Stigma-specific STIG1-like protein",
    str_detect(annotate_sp, "Biotin carboxyl carrier protein of acetyl-CoA carboxylase") ~ "Biotin carboxyl carrier protein of acetyl-CoA carboxylase",
    
    str_detect(annotate_sp, "heat shock protein") ~ "Heat shock protein domains",
    
    
    str_detect(annotate_sp, "Myb-related") |
      str_detect(annotate_sp, "MYB") |
      str_detect(annotate_sp, "MYB-like transcription factor") |
      str_detect(annotate_sp, "MYB31 transcription") ~ "MYB-related Transcription factors",

    str_detect(annotate_sp, "Riboflavin_synthase_like") |
      str_detect(annotate_sp, "ribE: riboflavin synthase, alpha subunit") |
      str_detect(annotate_sp, "Riboflavin_synthase_like") |
      str_detect(annotate_sp, "RIBOFLAVIN SYNTHASE") |
      # str_detect(annotate_sp, "RIBOFLAVIN SYNTHASE") |
      str_detect(annotate_sp, "Riboflavin synthase") ~ "Riboflavin synthase domain-like",

    str_detect(annotate_sp, "Universal stress protein")  ~ "Universal stress proteins",
    str_detect(annotate_sp, "ranscriptional regulator")  ~ "Transcriptional regulators (TAC1/EARS)",
    str_detect(annotate_sp, "Alpha/beta hydrolase domain-containing protein")  ~ "Alpha/beta hydrolase domain-containing protein",
    str_detect(annotate_sp, "Subtilisin-like protease")  ~ "Subtilisin-like proteases",
    str_detect(annotate_sp, "Dolichyl-phosphate beta-glucosyltransferase")  ~ "Dolichyl-phosphate beta-glucosyltransferases",
    # str_detect(annotate_sp, "Uncharacterized mitochondrial protein")  ~ "Uncharacterized mitochondrial proteins",
    # str_detect(annotate_sp, "Alpha/beta hydrolase domain-containing protein")  ~ "Alpha/beta hydrolase domain-containing protein",
    str_detect(annotate_sp, "OXIDOREDUCTASE, 2OG-FE II  OXYGENASE FAMILY PROTEIN") ~ "2OG-Fe(II) oxygenase superfamily",
    str_detect(annotate_sp, "DNA \\(cytosine-5\\)-methyltransferase \\w+")  ~ "DNA (cytosine-5)-methyltransferase (DRMB/DRM1/DRM2/DRM3)",
    str_detect(annotate_sp, "Sulfate/thiosulfate import ATP-binding protein")  ~ "Sulfate/thiosulfate import ATP-binding protein",
    str_detect(annotate_sp, "PHD finger protein") ~ "PHD finger proteins (MALE STERILITY 1/MEIOCYTE DEATH 1)",
    
    str_detect(annotate_sp, "\\(S\\)-8-oxocitronellyl enol synthase")  ~ "(S)-8-oxocitronellyl enol synthases",
    str_detect(annotate_sp, "Zinc finger protein") ~  "Zink finger proteins",
    str_detect(annotate_sp, "Uncharacterized") ~  "Uncharacterized protiens",
    str_detect(annotate_sp, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
    str_detect(annotate_sp, "PLASMODESMATA CALLOSE-BINDING PROTEIN") ~  "PLASMODESMATA CALLOSE-BINDING PROTEIN",
    str_detect(annotate_sp, "Transcription factor") ~  "Transcription factors (GTE7/GTE2)",
    str_detect(annotate_sp, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
    # str_detect(annotate_sp, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
    # str_detect(annotate_sp, "E3 ubiquitin-protein ligase") ~  "E3 ubiquitin-protein ligases",
    # 
    
    str_detect(annotate_sp, "Retrovirus-related") |
      str_detect(annotate_sp, "\\wransposon") |
      str_detect(annotate_sp, "Enzymatic polyprotein") |
      str_detect(annotate_sp, "Polyprotein") |
      str_detect(annotate_sp, "CHROMO DOMAIN-CONTAINING PROTEIN") |
      str_detect(annotate_sp, "MuDR family transposase") |
      str_detect(annotate_sp, "RNase_HI_RT_Ty3") |
      str_detect(annotate_sp, "scriptase") |
      str_detect(annotate_sp, "polymerases") |
      str_detect(annotate_sp, "COPIA") |
      str_detect(annotate_sp, "Genome polyprotein") |
      str_detect(annotate_sp, "DNA BINDING PROTEIN") |
      str_detect(annotate_sp, "POLYPROTEIN") |
      str_detect(annotate_sp, "CHROMO") |
      str_detect(annotate_sp, "Copia protein") |
      str_detect(annotate_sp, "transposable") |
      str_detect(annotate_sp, "Enzymatic polyprotein") |
      str_detect(annotate_sp, "NAC domain-containing protein") |
      
      str_detect(annotate_sp, "SPOSON") ~ "Transposable elements (TE)/TE domains",
    TRUE ~ annotate_sp),
  annotate_sp = str_remove(annotate_sp, ";")
  
)  %>%

dplyr::filter(start_int != "NA") %>%
  filter(
    signature_description != "consensus disorder prediction",
    signature_description != "Coil",
    signature_description != "OS06G0165300 PROTEIN",
    signature_description != "OS02G0307900 PROTEIN",
    signature_description != "OS01G0588500 PROTEIN-RELATED",
    
    signature_accession != "G3DSA:1.10.287.1490",
    signature_accession != "G3DSA:1.10.340.70",
    signature_accession != "G3DSA:1.20.5.340",
    signature_description != "Tropomyosin"
    # seq_md5 != "Reactome: R-HSA-8950505,"
  ) %>%
  distinct()  

EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 <- 
EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v1 %>%
# select(seq.name,stitle_bp_annot,stitle_sp,annotate_sp,annotate_intpr) %>%
  mutate(
    comb = case_when(str_detect(stitle_bp_annot,"hypothetical proteins") &
                       str_detect(stitle_sp,"NA")  ~ annotate_intpr,
                     # !str_detect(stitle_sp,"NA") &
                     !str_detect(annotate_sp,"NA")  ~ annotate_sp ),
    comb = str_replace_na(comb, "NA"),
    comb = case_when(comb=="NA" & str_detect(annotate_intpr,"\\w+") ~ annotate_intpr,TRUE ~ comb),
    comb = str_replace_na(comb, "NA"),
    comb = case_when(comb=="-" ~ "NA",TRUE ~ comb),
    comb = case_when(comb=="NA" & str_detect(stitle_bp_annot,"\\w+") ~ stitle_bp_annot,
                     TRUE ~ comb),
    
    comb = case_when(str_detect(comb, "SMALL NUCLEAR RNA ACTIVATING") | 
                       str_detect(comb,"SNRNA-ACTIVATING PROTEIN") |
                       str_detect(comb,"SMALL NUCLEAR RNA ACTIVATING COMPLEX") |
                       str_detect(comb,"Small nuclear RNA activating complex") ~ "Small nuclear RNA activating complex",
                     str_detect(comb, "E3 UBIQUITIN")  |
                       str_detect(comb, "UBIQUITIN-PROTEIN LIGASE") |
                       str_detect(comb, "E3 UBIQUITIN-PROTEIN") ~  "E3 ubiquitin proteins ligases",
                     str_detect(comb, "Mitochondrial carrier") ~  "Mitochondrial carrier protein",
                     str_detect(comb, "Jacalin-type lectin") ~  "Jacalin-related lectin",
                     str_detect(comb, "Chaperone protein DnaJ") ~  "Chaperone protein DnaJ",
                     str_detect(comb, "LEUCINE-RICH REPEAT (LRR) FAMILY PROTEIN") ~  "Leucine-rich repeat-containing protein",
                     str_detect(comb, "Pathogenesis-related thaumatin-like protein") | 
                       str_detect(comb, "Disease resistance proteins") |
                       str_detect(comb, "PR5-like receptor kinase") |
                       str_detect(comb, "PR5-like receptor kinase") |
                       str_detect(comb, "NDR1/HIN1-like protein") |
                       str_detect(comb, "Glucan endo-1,3-beta-glucosidase") |
                       str_detect(comb, "Pathogenesis-related") |
                       str_detect(comb, "Protein INCREASED RESISTANCE TO MYZUS PERSICAE") ~  "Disease resistance proteins (RPP8-like/RGAs/PRs)",
                     str_detect(comb, "\\wrotein of unknown function") |
                       str_detect(comb, "Domain of unknown function") |
                       str_detect(comb, "UPF0481 protein") ~  "Protein of unknown function",
                     str_detect(comb, "5'-NUCLEOTIDASE DOMAIN-CONTAINING") ~  "5' nucleotidase family",
                     # str_detect(comb, "5'-NUCLEOTIDASE DOMAIN-CONTAINING") ~  "5' nucleotidase family",
                     str_detect(comb, "26S prot") ~  "26S proteasome regulatory",
                     str_detect(comb, "Uncharacterized protiens") ~  "uncharacterized proteins",
                     str_detect(comb, "HMG boxes") |
                       str_detect(comb, "HMG-box domain") |
                       str_detect(comb, "HMG-box") ~ "High mobility group box domain",
                     str_detect(comb, "\\wucleoredoxin") ~  "Nucleoredoxin",
                     str_detect(comb, "NUCLEAR MIGRATION PROTEIN NUDC") ~  "NUCLEAR MOVEMENT PROTEIN NUDC",
                     str_detect(comb, "BSD domain") |
                       str_detect(comb, "BSD DOMAIN") ~  "BSD domain-like",
                     str_detect(comb, "ALUMINUM INDUCED PROTEIN") ~  "Aluminium induced protein",
                     str_detect(comb, "P-loop containing nucleoside triphosphate") ~  "P-loop containing nucleoside triphosphate hydrolase",
                     str_detect(comb, "TRANSLATION ELONGATION") ~  "Translational elongation factor related protein",
                     str_detect(comb, "ARM REPEAT SUPERFAMILY PROTEIN") |
                       str_detect(comb, "ARM repeat") ~  "ARM REPEAT SUPERFAMILY PROTEIN",
                     str_detect(comb, "CHROMO DOMAIN-CONTAINING PROTEIN") |
                       str_detect(comb, "Reverse transcriptase") |
                       str_detect(comb, "DNA-binding pseudobarrel domain") |
                       str_detect(comb, "Ribonuclease H") |
                       str_detect(comb, "retrotransposon") |
                       str_detect(comb, "reverse transcriptase") |
                       str_detect(comb, "retrotransposable") |
                       str_detect(comb, "Putative enzymatic polyprotein") |                      
                       str_detect(comb, "Integrase catalytic domain-containing protein") |
                       str_detect(comb, "RNA-DIRECTED DNA POLYMERASE HOMOLOG") ~  "Transposable elements (TE)/TE domains",
                     
                     str_detect(comb, "proline-rich") |
                       str_detect(comb, "basic proline-rich protein-like") |
                       str_detect(comb, " Proline rich extensin signature") ~  "Proline-rich protein",
                     str_detect(comb, "\\wlutamine-rich protein") ~  "Glutamine-rich protein 2-like",
                     
                     str_detect(comb, "Zink finger proteins") |
                       str_detect(comb, "Zinc finger CCHC-type superfamily") ~  "Zinc finger proteins",
                     # str_detect(comb, "ALUMINUM INDUCED PROTEIN") ~  "Aluminium induced protein",
                     # str_detect(comb, "ALUMINUM INDUCED PROTEIN") ~  "Aluminium induced protein",
                     
                     str_detect(comb, "DNA  CYTOSINE-5- -METHYLTRANSFERASE 3-RELATED") ~  "DNA (CYTOSINE-5)-METHYLTRANSFERASE DRM",
                     str_detect(comb, "Universal stress proteins") ~  "Universal stress protein family",
                     # str_detect(comb, "Universal stress proteins") ~  "Universal stress protein family",
                     str_detect(comb, "\\wentapeptide") ~  "Pentapeptide repeat-containing proteins",
                     
                     TRUE ~ comb)
    
    
    
    ) %>%
  
  mutate(
    comb = case_when(
      # str_detect(comb, "uncharacterized") |
      str_detect(comb, "unnamed protein product") |
        str_detect(comb, "uncharacterized protein") |
        str_detect(comb, "Protein of unknown function") |
        str_detect(comb, "DUF") |
        # str_detect(comb, "DUF") |
        str_detect(comb, "hypothetical proteins") ~  "Proteins of unknown function (unnamed protein product\n or uncharacterized proteins)",
      
      str_detect(comb, "Zink finger domains")  ~  "Zinc finger proteins",
      # str_detect(comb, "Pentapeptide repeat-like")  ~  "Pentapeptide repeat-containing proteins",
      # str_detect(comb, "\\entapeptide repeat-containing protein")  ~  "Pentapeptide repeat-containing proteins",
      str_detect(comb, "ARM REPEAT SUPERFAMILY PROTEIN")  ~  "ARM repeat superfamily protein",
      str_detect(comb, "MEIOSIS-SPECIFIC PROTEIN ASY2-RELATED")  ~  "Meiosis-specific protein ASY2-relateds",
      str_detect(comb, "DNA DAMAGE-INDUCIBLE PROTEIN 1-LIKE")  ~  "DNA damage-inducible protein 1-like",
      str_detect(comb, "GLYCINE-RICH CELL WALL STRUCTURAL PROTEIN 1.8-LIKE")  ~  "Glycine-rich cell wall structural protein 1.8-like",
      str_detect(comb, "MYOSIN-11-LIKE")  ~  "Myosin-11-like",
      str_detect(comb, "BAR")  ~  "BAR domain proteins",
      str_detect(comb, "AF4/FMR2 family member")  ~  "AF4/FMR2 family proteins",
      # str_detect(comb, "Disease resistance proteins \\(RPP8\\-like\\/RGAs\\/PRs\\)")  ~  "Disease resistance proteins (RPP8-like/RGAs/PRs)",
      
      TRUE ~ comb)
  ) %>%
  distinct()



# EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v3 <-
EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
select(seq.name,stitle_bp_annot,stitle_sp,annotate_sp,annotate_intpr,comb) %>%
  mutate(genome = case_when(str_detect(seq.name,"EVMZ") ~ "Mazia",
                            TRUE ~ "Bedadeti")) %>%
  # filter(genome != 'Mazia') %>%
  select(seq.name,genome,comb) %>%
  distinct() %>%
  group_by(genome,comb) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  arrange(desc(count)) %>% 
  as.data.frame() %>%
  filter(!str_detect(comb,"^OS"),
         !str_detect(comb,"PGCA protein"),
         !str_detect(comb,"Protein SON"),
         !str_detect(comb,"prisilkin-39-like"),
         !str_detect(comb,"Prokaryotic"),
         !str_detect(comb,"Glutamine-rich"),
         !str_detect(comb,"AF4/FMR2 family proteins"),
         # !str_detect(comb,"Prokaryotic"),
         # !str_detect(comb,"Prokaryotic"),
         # !str_detect(comb,"Prokaryotic"),
         # !str_detect(comb,"Prokaryotic"),
         # !str_detect(comb,"Prokaryotic"),
         # !str_detect(comb,"Prokaryotic"),
         count > 1 | 
           comb=="Heat shock protein domains" |
           comb=="Wound-induced basic protein "
           
        ) #%>%
  filter(str_detect(comb,"MATE"))




EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v3 %>%
  # select(seq.name,stitle_bp_annot,stitle_sp,annotate_sp,annotate_intpr,comb) %>%
  head()

# Mazia  
mazia_EV_specific_NOGO_plot  <-
EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v3 %>%
  dplyr::filter(genome == "Mazia") %>%
  ggplot(aes(x=reorder(comb, count), #Reorder gene sets by k/K values
             y= count)) +
  geom_col( width = 0.8) +
  # geom_text(aes(label = Count), hjust=1.1, vjust= 0.5, size=2, fontface = "bold")+
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  scale_y_continuous(breaks =c(2, 200,400,600,750))+
  # scale_x_discrete(position = "top") +
  coord_flip() +
  #fix labels
  labs(y="Proportion of significant genes over<br>the total genes" ,
       x="GO terms description for gene set (Specific to MABS specific)",
       #title = "Musa species specific" 
  ) +
  scale_fill_gradientn(colours = terrain.colors(10),
                       name  = "log10\npvalue")+
  
  # scale_fill_continuous(name  = "log10 pvalue")+
  theme( 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 6,face = "bold"),
    axis.text.x = element_text(size = 6),
    # axis.title.x
    # panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    # panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_text(size=6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.35,"cm"),
    legend.position = "top"
    # legend.position = c(0.7,0.85)
  )

## Bedadeti 

bedadeti_EV_specific_NOGO_plot  <-
  EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v3 %>%
  dplyr::filter(genome != "Mazia") %>%
  # dplyr::filter(frac_cov==0.30) %>%
  # dplyr::filter(p.adjust < 0.05) %>%
  # filter(genome=="*E. ventricosum*") %>%
  # dplyr::filter(description != "NA") %>% 
  # dplyr::filter(aspects == 'Biological process') %>%
  # dplyr::filter(p.adjust < 0.05) %>%
  # filter(Significant <= 2) %>%
  # filter(!str_detect(description, "^protein transport to"), 
  #        !str_detect(description, "^positive regulation"), 
  #        !str_detect(description, "^peptidyl")) %>%
  # filter(Category == "EV shared") %>% 
ggplot(aes(x=reorder(comb, count), #Reorder gene sets by k/K values
           y= count)) +
  geom_col( width = 0.8) +
  # geom_text(aes(label = Count), hjust=1.1, vjust= 0.5, size=2, fontface = "bold")+
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  scale_y_continuous(breaks =c(2,20,40,58))+
  scale_x_discrete(position = "top") +
  coord_flip() +
  #fix labels
  labs(y="Proportion of significant genes over<br>the total genes" ,
       x="GO terms description for gene set (Specific to MABS specific)",
       #title = "Musa species specific" 
  ) +
  scale_fill_gradientn(colours = terrain.colors(10),
                       name  = "log10\npvalue")+
  
  # scale_fill_continuous(name  = "log10 pvalue")+
  theme( 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_markdown(face = "bold"),
    axis.text.y = element_text(size = 6,face = "bold"),
    axis.text.x = element_text(size = 6),
    # axis.title.x
    # panel.grid.major.x  = element_line(color ="#888888", size = 0.08),
    # panel.background = element_rect(fill="#FFFFFF", color=NA),
    legend.title = element_text(size=6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.35,"cm"),
    legend.position = "top"
    # legend.position = c(0.7,0.85)
  )
# library(patchwork)
(mazia_EV_specific_NOGO_plot|bedadeti_EV_specific_NOGO_plot)
ggsave("../go_terms_enrichment_v2/EV_mazia_bedadeti_specific.NOGO.jpeg", width=7.5, height=4)



Ensete_MUSA_seq %>%
  left_join(
EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
  select(seq.name,seq.text,comb) %>%
  bind_rows(MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro.v2  %>%
  select(seq.name,seq.text,comb))) %>%
  mutate(comb = str_replace_na(comb,"NA")) %>%
  filter(!str_detect(comb,"^Transposable elements"),
         !str_detect(comb,"^Proteins of unknown function")) %>%
  select(seq.name,seq.text) %>%
  distinct() %>% 
  mutate(seq.name = paste0(">",seq.name)) %>%
  pivot_longer(cols = c(seq.name,seq.text)) %>%
  select(value) %>%
  write.table(paste0(panev_path_out_v1, "EV_MZBD_MAB.prot.non_TE.fa"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)


## MA CDS 

Ensete_MUSA_seq %>%
  left_join(
    EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
      select(seq.name,seq.text,comb) %>%
      bind_rows(MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro.v2  %>%
                  select(seq.name,seq.text,comb))) %>%
  mutate(comb = str_replace_na(comb,"NA")) %>%
  filter(!str_detect(comb,"^Transposable elements"),
         !str_detect(comb,"^Proteins of unknown function")) %>%
  select(seq.name,seq.text) %>%
  filter(!str_detect(seq.name,"EV")) %>%
  left_join(Musa_AA_BB_gene_CDS_gene) %>%
  filter(str_detect(ref_genome,"MA")) %>%
  mutate(len = end - start) %>%
  filter(len > 160 ) %>% 
  select(-len) %>%
  distinct() %>%
  filter(str_detect(contig_chr_name,"^\\w+")) %>%
  select(contig_chr_name,start,end) %>%  
  write.table(paste0(panev_path_out_v1, "MA.non_TE.CDS.bed"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)



## MB CDS 

Ensete_MUSA_seq %>%
  left_join(
    EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
      select(seq.name,seq.text,comb) %>%
      bind_rows(MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro.v2  %>%
                  select(seq.name,seq.text,comb))) %>%
  mutate(comb = str_replace_na(comb,"NA")) %>%
  filter(!str_detect(comb,"^Transposable elements"),
         !str_detect(comb,"^Proteins of unknown function")) %>%
  select(seq.name,seq.text) %>%
  filter(!str_detect(seq.name,"EV")) %>%
  left_join(Musa_AA_BB_gene_CDS_gene) %>%
  filter(str_detect(ref_genome,"MB")) %>%
  select(contig_chr_name,start,end) %>%
  mutate(len = end - start) %>%
  filter(len > 160 ) %>% 
  select(-len) %>%
  distinct() %>%
  filter(str_detect(contig_chr_name,"^\\w+")) %>%
  select(contig_chr_name,start,end) %>%
  write.table(paste0(panev_path_out_v1, "MB.non_TE.CDS.bed"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)



## EV mazia 

Ensete_MUSA_seq %>%
  left_join(
    EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
      select(seq.name,seq.text,comb) %>%
      bind_rows(MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro.v2  %>%
                  select(seq.name,seq.text,comb))) %>%
  mutate(comb = str_replace_na(comb,"NA")) %>%
  filter(!str_detect(comb,"^Transposable elements"),
         !str_detect(comb,"^Proteins of unknown function")) %>%
  select(seq.name,seq.text) %>%
  filter(str_detect(seq.name,"EV")) %>%
  select(-seq.text) %>%
  # head() 
  left_join(
    mazia_bedadeti_gff %>%
      filter(type == "CDS") %>%
      mutate(seq.name = str_extract(attributes,"Parent=\\w+.*[^;]"),
             seq.name = str_remove(seq.name,"Parent=")) %>%
      left_join(EVMZBD_full_length_CDS_gene_ID) %>%
      select(seqid,start,end,seq.name)
  ) %>%
  filter(str_detect(seq.name,"EVMZ")) %>%
  select(seqid,start,end) %>%
  mutate(len = end - start) %>%
  filter(len > 160 ) %>% 
  select(-len) %>%
  filter(str_detect(seqid,"^\\w+")) %>%
  distinct() %>%
  write.table(paste0(panev_path_out_v1, "EVMZ.non_TE.CDS.bed"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)



## EVBD


Ensete_MUSA_seq %>%
  left_join(
    EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro.v2 %>%
      select(seq.name,seq.text,comb) %>%
      bind_rows(MAB_specific_blasp_ncbi_uniprot_25_perc_CDS_max_intpro.v2  %>%
                  select(seq.name,seq.text,comb))) %>%
  mutate(comb = str_replace_na(comb,"NA")) %>%
  filter(!str_detect(comb,"^Transposable elements"),
         !str_detect(comb,"^Proteins of unknown function")) %>%
  select(seq.name,seq.text) %>%
  filter(str_detect(seq.name,"EV")) %>%
  select(-seq.text) %>%
  # head() 
  left_join(
    mazia_bedadeti_gff %>%
      filter(type == "CDS") %>%
      mutate(seq.name = str_extract(attributes,"Parent=\\w+.*[^;]"),
             seq.name = str_remove(seq.name,"Parent=")) %>%
      left_join(EVMZBD_full_length_CDS_gene_ID) %>%
      select(seqid,start,end,seq.name)
  ) %>%
  filter(!str_detect(seq.name,"EVMZ")) %>%
  mutate(len = end - start) %>%
  filter(len > 160 ) %>% 
  select(-len) %>%
  filter(str_detect(seqid,"^\\w+")) %>%
  distinct() %>%
  filter(str_detect(seqid,"^\\w+")) %>%
  select(seqid,start,end) %>%
  # head()
  write.table(paste0(panev_path_out_v1, "EVBD.non_TE.CDS.bed"), # name of the fasta file to be saved
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              na = "",
              quote=FALSE)

  
  
  library(ape)
  
  mazia_bedadeti_gff <- 
  read.gff("~/rstudio/bgimazia_musa/reannotation_analysis/manuscript/Supplementary_files/bedadeti_maker_annotation.gff3") %>%
    rbind(read.gff("~/rstudio/bgimazia_musa/reannotation_analysis/manuscript/Supplementary_files/mazia_maker_annotation.gff3"))
    
  mazia_bedadeti_gff %>%
    filter(type == "CDS") %>%
    mutate(seq.name = str_extract(attributes,"Parent=\\w+.*[^;]"),
           seq.name = str_remove(seq.name,"Parent=")) %>%
    left_join(EVMZBD_full_length_CDS_gene_ID) %>%
    select(seqid,start,end,seq.name)
    
  head()

  Musa_AA_BB_gene_CDS_gene %>%
  tail()

Musa_AA_BB_gene_CDS

EVMZBD_full_length_CDS_gene_ID %>%
  head()














# head()

#%>%
  dplyr::filter(str_detect(annotate_sp, "(S)-8-oxocitronellyl enol synthase"))
  
  head()
  head()
  colnames()
  # # arrange() 

  # mutate(
  #   
  #   annotate_intpr = case_when(str_detect(annotate_intpr,"PROLINE-RICH EXTENSIN") ~  "Proline rich extensin signature" )
  # )
  # 
  
  write.table(paste0(panev_path_out_v1,"EV_specifi_interproscan_NOGO_blastp_stitile.txt"), 
              col.names = T,
              row.names = F,
              sep = '\t',
              quote = F)




filePath = paste0(panev_path_out_v1,"EV_specifi_interproscan_NOGO_blastp_stitile.txt")
# text <- readLines(filePath)

text_EV <-
  read.delim(filePath, header = T) %>%
  mutate(stitle = str_replace(stitle,"Ensete ventricosum hypothetical proteins", "EV hypothetical proteins")) %>%
  mutate(stitle = str_replace(stitle, "RM repeat protein interacting WITH ABF2 protein", "RM repeat protein interacting WITH ABF2")) %>%
  # mutate( stitle = str_to_upper(stitle)) %>%
  group_by(stitle) %>%
  summarise(frq = n()) %>%
  ungroup() %>%
  arrange(desc(frq))


text %>%
  head()
wordcloud(words  = text_EV$stitle, 
          freq = text_EV$frq, 
          min.freq = 1,
          max.words=200, 
          random.order=FALSE, 
          rot.per=0.1, 
          colors=brewer.pal(8, "Dark2"))

ggsave()

ggsave("PanMusa_wordcloud.jpeg", width=5.6, height=4)


png("PanMusa_wordcloud.v2.png", width=12,height=8, units='in', res=300)
par(mar = rep(0, 4))
wordcloud(names(tb),as.numeric(tb), 
          # scale=c(8,.3),
          min.freq=1,
          max.words=100, 
          andom.order=T, 
          rot.per=0.2, 
          random.order=FALSE,
          olors=brewer.pal(8, "Dark2"))



