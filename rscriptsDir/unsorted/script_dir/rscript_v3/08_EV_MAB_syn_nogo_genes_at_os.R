EV_MAB_nogo_noreads_at_os_nr <-
read.delim("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir//EV_MAB.nogo_noreads.atha_osa.diamond.blastp.e05.out.gz") %>% 
  # EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro %>%
  # mutate(stitle_sp = str_replace_na(stitle_sp, "NA")) %>%
  distinct() %>%
  mutate(
    stitle_bp_annot = stitle,
    stitle_bp_annot = str_replace(stitle_bp_annot, "isoform \\w+\\d+", "isoform"),
    # stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+\\s+\\|\\s+\\w+:\\s+"),
    # stitle_bp_annot = str_remove(stitle_bp_annot,"\\s+\\w+=\\w+\\s+.*"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+.\\d+\\s+\\|\\s+\\w+:\\s+\\w+\\w+\\w+\\w+"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"\\s+\\|\\schr\\d+:\\d+-\\d+\\s+\\w+.*"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"LOC_\\w+.\\d+\\s+\\w+\\|"),
    stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+.\\d+\\s+\\|\\s+\\w+:\\s+"),
    stitle_bp_annot = str_remove(stitle_bp_annot,",\\s+\\w+,\\s+\\w+,\\s+\\w+"),
    stitle_bp_annot = str_remove(stitle_bp_annot,".*\\|"),
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
      str_detect(stitle_bp_annot, "heat shock protein") |
        str_detect(stitle_bp_annot, "HEAT SHOCK") ~ "Heat shock protein domains",
      
      str_detect(stitle_bp_annot, "Myb-related") |
        str_detect(stitle_bp_annot, "MYB") |
        str_detect(stitle_bp_annot, "MYB-like transcription factor") |
        str_detect(stitle_bp_annot, "myb domain protein") |
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
      str_detect(stitle_bp_annot, "\\w+inc-finger homeodomain.*") |
        str_detect(stitle_bp_annot, "ZINC FINGER HOMEODOMAIN") ~ "Zinc-finger homeodomain protein",
      str_detect(stitle_bp_annot, "DNA-directed RNA polymera")  ~ "DNA-directed RNA polymerase II",
      str_detect(stitle_bp_annot, "RNA-directed DNA polymera")  ~ "RNA-directed DNA polymerase",
      str_detect(stitle_bp_annot, "keratin-associated protein")  ~ "keratin-associated proteins",
      str_detect(stitle_bp_annot, "mRNA turnover protein.*")  ~ "mRNA turnover proteins",
      str_detect(stitle_bp_annot, "zf-CCHC domain-containing protein.*")  ~ "zf-CCHC domain-containing proteins",
      str_detect(stitle_bp_annot, "wall-associated receptor kinase.*")  ~ "wall-associated receptor kinases",
      str_detect(stitle_bp_annot, "RPM1-interacting")  ~ "RPM1-interacting protein 4",
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
      str_detect(stitle_bp_annot, "60S ribosomal.*")  |
        str_detect(stitle_bp_annot, "60S acidic ribosomal.*") ~  "60S ribosomal proteins",
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
        str_detect(stitle_bp_annot, "\\w+acalin") |
        str_detect(stitle_bp_annot, "\\w+acalin-related") ~ "Jacalin-related lectin",
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
        str_detect(stitle_bp_annot, "RNA-DIRECTED DNA POLYMERASE HOMOLOG") ~  "Transposable elements (TE)/TE domains" ,
      
      str_detect(stitle_bp_annot, "proline-rich") |
        str_detect(stitle_bp_annot, "basic proline-rich protein-like") |
        str_detect(stitle_bp_annot, " Proline rich extensin signature") ~  "Proline-rich protein",
      str_detect(stitle_bp_annot, "\\w+lutamine-rich protein") ~  "Glutamine-rich protein 2-like",
      
      str_detect(stitle_bp_annot, "Zink finger proteins") |
        str_detect(stitle_bp_annot, "Zinc finger CCHC-type superfamily") ~  "Zinc finger proteins",
      str_detect(stitle_bp_annot, "DNA-damage-repair.*") |
        str_detect(stitle_bp_annot, "DNA damage-repair.*")  | 
        str_detect(stitle_bp_annot, "DNA-DAMAGE REPAIR")~  "DNA-damage-repair/toleration protein",
      
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
      str_detect(stitle_bp_annot, "chain acyl-CoA synthetase") |
        str_detect(stitle_bp_annot, "LONG-CHAIN ACYL-COA SYNTHASE") |
        str_detect(stitle_bp_annot, "Long-Chain Acyl-Coenzyme") ~  "Long chain acyl-CoA synthetase",
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
      str_detect(stitle_bp_annot, "ACETYLTRANSFERASE-RELATED PROTEIN") ~  "P300/CBP ACETYLTRANSFERASE-RELATED PROTEIN",
      str_detect(stitle_bp_annot, "CAMK includes calcium") ~  "CAMK includes calcium/calmodulin depedent protein kinases",
      str_detect(stitle_bp_annot, "Flowering Locus T") ~  "FT-Like5 homologous to Flowering Locus T  gene",
      str_detect(stitle_bp_annot, "phospholipase") ~  "phospholipase D",
      str_detect(stitle_bp_annot, "Levadura ") ~  "Arabidopsis TÃ³xicos en Levadura",
      str_detect(stitle_bp_annot, "hosphatidylethanolamine") ~  "hosphatidylethanolamine-binding protein",
      str_detect(stitle_bp_annot, "ATP-binding cassette") ~  "ATP-binding cassette proteins",
      str_detect(stitle_bp_annot, "DNA J protein") ~  "DNA J protein proteins",
      str_detect(stitle_bp_annot, "attachment domain containing") ~  "biotin/lipoyl attachment domain containing",
      str_detect(stitle_bp_annot, "gamma carbonic anhydrase") ~  "gamma carbonic anhydrase",
      str_detect(stitle_bp_annot, "NDR") ~  "NDR or NDR1/HIN1 like",
      str_detect(stitle_bp_annot, "\\w+erpene synthase ") ~  "terpene synthase",
      str_detect(stitle_bp_annot, "Hydrolase Domain-containing Protein") ~  "ABHD17 (Alpha/Beta Hydrolase Domain-containing Protein 17)-like Acyl Protein Thioesterase",
      str_detect(stitle_bp_annot, "ACT domain repeats") ~  "ACT domain repeats",
      str_detect(stitle_bp_annot, "Arabidopsis toxicos en levadura") ~  "Arabidopsis toxicos en levadura",
      str_detect(stitle_bp_annot, "EXORDIUM like") ~  "EXORDIUM like",
      str_detect(stitle_bp_annot, "FORKED-LIKE") ~  "FORKED-LIKE",
      str_detect(stitle_bp_annot, "regulatory components of ABA receptor") ~  "regulatory components of ABA receptor",
      str_detect(stitle_bp_annot, "hypersensitive induced reaction") ~  "hypersensitive induced reaction",
      str_detect(stitle_bp_annot, "ATP-dependent Clp protease ATP-binding subunit") ~  "ATP-dependent Clp protease ATP-binding subunit",
      str_detect(stitle_bp_annot, "RECOGNITION OF PERONOSPORA PARASITICA") ~  "RECOGNITION OF PERONOSPORA PARASITICA",
      str_detect(stitle_bp_annot, "chrc:110398-112638") ~  "Chloroplast encoded NADH dehydrogenase unit",
      str_detect(stitle_bp_annot, "chrc:152806-154312") ~  "RPL2.1 encodes a chloroplast ribosomal protein L2",
      str_detect(stitle_bp_annot, "chrc:84337-8584") ~  "RPL2.1 encodes a chloroplast ribosomal protein L2",
      str_detect(stitle_bp_annot, "chrc:9938-11461") ~  "ATP synthase subunit alpha and located in chloroplast",
      
      
      TRUE ~ stitle_bp_annot),
    
    
    # stitle_bp_annot = str_remove(stitle_bp_annot, "\\w+|\\w+|\\w+\\s+"),
    # stitle_bp_annot = str_remove(stitle_bp_annot, "PREDICTED:\\s+"),

    stitle_bp_annot = str_remove(stitle_bp_annot, ", putative, expressed"),
    
    
    
  ) %>% 
  # filter(str_detect(stitle_bp_annot,"Sex determination")) %>%
  mutate(cov = round(length/qlen,2)) %>%
  filter(cov >= 0.8) %>%
  filter(pident >= 0.5) %>%
  # distinct(stitle_bp_annot) %>%
  group_by(qseqid) %>%
  summarise(bitscore = max(bitscore)) %>%
  ungroup() %>%
  mutate(cat = "selected") %>%
  left_join(
    # EV_MAB_nogo_noreads_uniprot_nr <-
    read.delim("D:/persdoc_dir/allProject_sciptDir/local_scripts/ev_mab_comparative_genomics/data_dir//EV_MAB.nogo_noreads.atha_osa.diamond.blastp.e05.out.gz") %>% 
      # EV_specific_blasp_ncbi_uniprot_25_perc_CDS_max_interpro %>%
      # mutate(stitle_sp = str_replace_na(stitle_sp, "NA")) %>%
      distinct() %>%
      mutate(
        stitle_bp_annot = stitle,
        stitle_bp_annot = str_replace(stitle_bp_annot, "isoform \\w+\\d+", "isoform"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+\\s+\\|\\s+\\w+:\\s+"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"\\s+\\w+=\\w+\\s+.*"),
        stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+.\\d+\\s+\\|\\s+\\w+:\\s+\\w+\\w+\\w+\\w+"),
        stitle_bp_annot = str_remove(stitle_bp_annot,"\\s+\\|\\schr\\d+:\\d+-\\d+\\s+\\w+.*"),
        stitle_bp_annot = str_remove(stitle_bp_annot,"LOC_\\w+.\\d+\\s+\\w+\\|"),
        stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+.\\d+\\s+\\|\\s+\\w+:\\s+"),
        stitle_bp_annot = str_remove(stitle_bp_annot,",\\s+\\w+,\\s+\\w+,\\s+\\w+"),
        stitle_bp_annot = str_remove(stitle_bp_annot,".*\\|"),
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
          str_detect(stitle_bp_annot, "heat shock protein") |
            str_detect(stitle_bp_annot, "HEAT SHOCK") ~ "Heat shock protein domains",
          
          str_detect(stitle_bp_annot, "Myb-related") |
            str_detect(stitle_bp_annot, "MYB") |
            str_detect(stitle_bp_annot, "MYB-like transcription factor") |
            str_detect(stitle_bp_annot, "myb domain protein") |
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
          str_detect(stitle_bp_annot, "\\w+inc-finger homeodomain.*") |
            str_detect(stitle_bp_annot, "ZINC FINGER HOMEODOMAIN") ~ "Zinc-finger homeodomain protein",
          str_detect(stitle_bp_annot, "DNA-directed RNA polymera")  ~ "DNA-directed RNA polymerase II",
          str_detect(stitle_bp_annot, "RNA-directed DNA polymera")  ~ "RNA-directed DNA polymerase",
          str_detect(stitle_bp_annot, "keratin-associated protein")  ~ "keratin-associated proteins",
          str_detect(stitle_bp_annot, "mRNA turnover protein.*")  ~ "mRNA turnover proteins",
          str_detect(stitle_bp_annot, "zf-CCHC domain-containing protein.*")  ~ "zf-CCHC domain-containing proteins",
          str_detect(stitle_bp_annot, "wall-associated receptor kinase.*")  ~ "wall-associated receptor kinases",
          str_detect(stitle_bp_annot, "RPM1-interacting")  ~ "RPM1-interacting protein 4",
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
          str_detect(stitle_bp_annot, "60S ribosomal.*")  |
            str_detect(stitle_bp_annot, "60S acidic ribosomal.*") ~  "60S ribosomal proteins",
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
            str_detect(stitle_bp_annot, "\\w+acalin") |
            str_detect(stitle_bp_annot, "\\w+acalin-related") ~ "Jacalin-related lectin",
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
            str_detect(stitle_bp_annot, "RNA-DIRECTED DNA POLYMERASE HOMOLOG") ~  "Transposable elements (TE)/TE domains" ,
          
          str_detect(stitle_bp_annot, "proline-rich") |
            str_detect(stitle_bp_annot, "basic proline-rich protein-like") |
            str_detect(stitle_bp_annot, " Proline rich extensin signature") ~  "Proline-rich protein",
          str_detect(stitle_bp_annot, "\\w+lutamine-rich protein") ~  "Glutamine-rich protein 2-like",
          
          str_detect(stitle_bp_annot, "Zink finger proteins") |
            str_detect(stitle_bp_annot, "Zinc finger CCHC-type superfamily") ~  "Zinc finger proteins",
          str_detect(stitle_bp_annot, "DNA-damage-repair.*") |
            str_detect(stitle_bp_annot, "DNA damage-repair.*")  | 
            str_detect(stitle_bp_annot, "DNA-DAMAGE REPAIR")~  "DNA-damage-repair/toleration protein",
          
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
          str_detect(stitle_bp_annot, "chain acyl-CoA synthetase") |
            str_detect(stitle_bp_annot, "LONG-CHAIN ACYL-COA SYNTHASE") |
            str_detect(stitle_bp_annot, "Long-Chain Acyl-Coenzyme") ~  "Long chain acyl-CoA synthetase",
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
          str_detect(stitle_bp_annot, "ACETYLTRANSFERASE-RELATED PROTEIN") ~  "P300/CBP ACETYLTRANSFERASE-RELATED PROTEIN",
          str_detect(stitle_bp_annot, "CAMK includes calcium") ~  "CAMK includes calcium/calmodulin depedent protein kinases",
          str_detect(stitle_bp_annot, "Flowering Locus T") ~  "FT-Like5 homologous to Flowering Locus T  gene",
          str_detect(stitle_bp_annot, "phospholipase") ~  "phospholipase D",
          str_detect(stitle_bp_annot, "Levadura ") ~  "Arabidopsis TÃ³xicos en Levadura",
          str_detect(stitle_bp_annot, "hosphatidylethanolamine") ~  "hosphatidylethanolamine-binding protein",
          str_detect(stitle_bp_annot, "ATP-binding cassette") ~  "ATP-binding cassette proteins",
          str_detect(stitle_bp_annot, "DNA J protein") ~  "DNA J protein proteins",
          str_detect(stitle_bp_annot, "attachment domain containing") ~  "biotin/lipoyl attachment domain containing",
          str_detect(stitle_bp_annot, "gamma carbonic anhydrase") ~  "gamma carbonic anhydrase",
          str_detect(stitle_bp_annot, "NDR") ~  "NDR or NDR1/HIN1 like",
          str_detect(stitle_bp_annot, "\\w+erpene synthase ") ~  "terpene synthase",
          str_detect(stitle_bp_annot, "Hydrolase Domain-containing Protein") ~  "ABHD17 (Alpha/Beta Hydrolase Domain-containing Protein 17)-like Acyl Protein Thioesterase",
          str_detect(stitle_bp_annot, "ACT domain repeats") ~  "ACT domain repeats",
          str_detect(stitle_bp_annot, "Arabidopsis toxicos en levadura") ~  "Arabidopsis toxicos en levadura",
          str_detect(stitle_bp_annot, "EXORDIUM like") ~  "EXORDIUM like",
          str_detect(stitle_bp_annot, "FORKED-LIKE") ~  "FORKED-LIKE",
          str_detect(stitle_bp_annot, "regulatory components of ABA receptor") ~  "regulatory components of ABA receptor",
          str_detect(stitle_bp_annot, "hypersensitive induced reaction") ~  "hypersensitive induced reaction",
          str_detect(stitle_bp_annot, "ATP-dependent Clp protease ATP-binding subunit") ~  "ATP-dependent Clp protease ATP-binding subunit",
          str_detect(stitle_bp_annot, "RECOGNITION OF PERONOSPORA PARASITICA") ~  "RECOGNITION OF PERONOSPORA PARASITICA",
          str_detect(stitle_bp_annot, "chrc:110398-112638") ~  "Chloroplast encoded NADH dehydrogenase unit",
          str_detect(stitle_bp_annot, "chrc:152806-154312") ~  "RPL2.1 encodes a chloroplast ribosomal protein L2",
          str_detect(stitle_bp_annot, "chrc:84337-8584") ~  "RPL2.1 encodes a chloroplast ribosomal protein L2",
          str_detect(stitle_bp_annot, "chrc:9938-11461") ~  "ATP synthase subunit alpha and located in chloroplast",
          
          
          TRUE ~ stitle_bp_annot),
        
        
        stitle_bp_annot = str_remove(stitle_bp_annot, ", putative, expressed"),
        # stitle_bp_annot = str_remove(stitle_bp_annot, "PREDICTED:\\s+"),
        # 
        # stitle_bp_annot = str_remove(stitle_bp_annot,"^\\w+.\\d+\\s+"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"isoform.*"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"LOC\\d.*+"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+_\\d+.*"),
        # stitle_bp_annot = str_remove(stitle_bp_annot,"\\w+_\\w+.*"),
        
        
      ) %>% 
      # filter(str_detect(stitle_bp_annot,"Sex determination")) %>%
      mutate(cov = round(length/qlen,2)) %>%
      filter(cov >= 0.8) %>%
      filter(pident >= 0.5)
  ) %>%
mutate(cat = str_replace_na(cat,"NA")) %>%
  filter(cat != "NA") #%>%



## atos
EV_MAB_nogo_noreads_at_os_nr %>% 
  select(qseqid, stitle_bp_annot) %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) %>%
  # filter(str_detect(qseqid,"EVBD")) %>%
  filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
         stitle_bp_annot != "Uncharacterized proteins" ,
         stitle_bp_annot != "unnamed protein product",
         stitle_bp_annot != " no full name available",
         stitle_bp_annot != "unnamed protein product, partial") %>%
  group_by(stitle_bp_annot) %>%
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(desc(count))

## ncbi nr

EV_MAB_nogo_noreads_ncbi_nr %>%
  select(qseqid, stitle_bp_annot) %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) %>%
  # filter(str_detect(qseqid,"EVMZ")) %>%
  filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
         stitle_bp_annot != "Uncharacterized proteins" ,
         stitle_bp_annot != "Unnamed protein product",
         stitle_bp_annot != " no full name available",
         stitle_bp_annot != "unnamed protein product, partial") %>%
  group_by(stitle_bp_annot) %>%
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(desc(count))

## uniprot

EV_MAB_nogo_noreads_uniprot_nr %>% 
  select(qseqid, stitle_bp_annot) %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) %>%
  filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
         stitle_bp_annot != "Uncharacterized proteins" ,
         stitle_bp_annot != "Unnamed protein product",
         stitle_bp_annot != " no full name available",
         stitle_bp_annot != "unnamed protein product, partial") %>%
  group_by(stitle_bp_annot) %>%
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(desc(count))




## Genes with homology in ncbi_nr which have also hits in uniprot 

EV_MAB_nogo_noreads_ncbi_uniprot_hits <-
  EV_MAB_nogo_noreads_ncbi_nr %>%
  select(qseqid, stitle_bp_annot) %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) %>%
  # filter(str_detect(qseqid,"EVMZ")) %>%
  filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
         stitle_bp_annot != "Uncharacterized proteins" ,
         stitle_bp_annot != "Unnamed protein product",
         stitle_bp_annot != " no full name available",
         stitle_bp_annot != "unnamed protein product, partial") %>%
  distinct() %>% 
  rename(stitle_bp_annot_ncbinr=stitle_bp_annot) %>%
  select(-genome)%>%
  as.data.frame() %>%
    full_join(
    # Find and join their mapping if these protein have blast hit match in uniprot curated proteins 
    # left_join(
    EV_MAB_nogo_noreads_uniprot_nr %>% 
      select(qseqid, stitle_bp_annot) %>%
      mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                                str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                                str_detect(qseqid,"Mac") ~ "MA",
                                str_detect(qseqid,"Mba") ~ "MB")) %>%
      filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
             stitle_bp_annot != "Uncharacterized proteins" ,
             stitle_bp_annot != "Unnamed protein product",
             stitle_bp_annot != " no full name available",
             stitle_bp_annot != "unnamed protein product, partial") %>%
      # distinct(qseqid) %>%
      as.data.frame() %>%
      mutate(db = "uniprot")) %>% 
    mutate(stitle_bp_annot_ncbinr = str_replace_na(stitle_bp_annot_ncbinr,"NA"),
           stitle_bp_annot = str_replace_na(stitle_bp_annot,"NA"),
           stitle_bp_annot_ncbinr = case_when(stitle_bp_annot_ncbinr=="NA" ~ stitle_bp_annot,
                                       TRUE ~ stitle_bp_annot_ncbinr)) %>%
    select(-stitle_bp_annot,-genome,-db) %>%
  rename(bp_ncbinr_uniprot = stitle_bp_annot_ncbinr)


  
## all combined blast hits of EVMAB genes that lack GO-terms annotation 

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits <- 
EV_MAB_nogo_noreads_at_os_nr %>% 
  select(qseqid, stitle_bp_annot) %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) %>%
  # filter(str_detect(qseqid,"EVBD")) %>%
  filter(stitle_bp_annot  !=  "Hypothetical proteins" ,
         stitle_bp_annot != "Uncharacterized proteins" ,
         stitle_bp_annot != "unnamed protein product",
         stitle_bp_annot != " no full name available",
         stitle_bp_annot != "unnamed protein product, partial") %>%
  full_join(EV_MAB_nogo_noreads_ncbi_uniprot_hits) %>%
  mutate(stitle_bp_annot = str_replace(stitle_bp_annot, "AMP-binding enzyme","Long chain acyl-CoA synthetase (AMP-binding enzyme)"),
         stitle_bp_annot = str_replace(stitle_bp_annot, "CAMK includes calcium/calmodulin depedent protein kinases","CAMK includes calcium/calmodulin depedent protein kinasese (MAPK)"),
         stitle_bp_annot = str_replace(stitle_bp_annot, "starch synthase","starch synthase (Chloroplastic proteins)")) %>%
  mutate(stitle_bp_annot = str_replace_na(stitle_bp_annot,"NA"),
         bp_ncbinr_uniprot = str_replace_na(bp_ncbinr_uniprot,"NA"),
         stitle_bp_annot = case_when(stitle_bp_annot=="NA" ~ bp_ncbinr_uniprot,
                                     TRUE ~ stitle_bp_annot)
         ) %>%
  select(-bp_ncbinr_uniprot,-genome) %>%
  distinct() %>%
  mutate(genome = case_when(str_detect(qseqid,"EVMZ") ~ "EV_mazia",
                            str_detect(qseqid,"EVBD") ~ "EV_bedadeti",
                            str_detect(qseqid,"Mac") ~ "MA",
                            str_detect(qseqid,"Mba") ~ "MB")) 



  

# double check that all of the genes with blast hits do not have go terms 

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits <-
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  as.data.frame() %>%
  separate(qseqid,into = c("gene_id","isomer"), sep = "-") %>%
  # filter(str_detect(gene_id,"^Ma")) %>%
  left_join(GO_terms_global_genes) %>%
  mutate(GO_terms = str_replace_na(GO_terms,"NA")) %>%
  filter(GO_terms == "NA") %>%
  select(-GO_terms) %>%
  distinct()



# extract EVMAB proteins/genes that lack GO-term annotation and showed either no blastp hists against public protein database or resulted 
# homology with Hypothetical/Uncharacterized/unmaed proteins

GO_terms_global_genes %>%
  # filter(str_detect(gene_id,"^Ma")) %>%
  head()

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
  head()

EV_MAB_specific_genes_nogo_noblastp <-
EV_MAB_custered_genes %>%
  filter(clust=="EV_bedadeti_specific" |
           clust=="EV_mazia_specific" |
           clust=="MA_specific" |
           clust=="MB_specific") %>%
  # join with genes who have go terms
  full_join(GO_terms_global_genes) %>%
  mutate(GO_terms = str_replace_na(GO_terms,"NA")) %>%
  filter(GO_terms == "NA") %>%
  select(-GO_terms) %>%
  # Filter those with blast hits
  full_join(EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits) %>%
  mutate(stitle_bp_annot  = str_replace_na(stitle_bp_annot, "NA")) %>%
  filter(stitle_bp_annot  == "NA") %>%
  select(gene_id,clust) %>%
  distinct() %>%
  rename(Name = gene_id) %>%
  # head() %>%
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      mutate(Name = case_when(str_detect(Name,"Ma") |
                                str_detect(Name,"Mb")~ str_c(Name,".1"),
                              TRUE ~ Name))
  ) %>%
  as.data.frame() %>%
  # mutate(isomers = str_replace_na(isomers, "")) %>%
  # mutate(Name = case_when(str_detect(Name,"-") ~ str_c(Name,isomers,sep = "-"),
  #                         TRUE ~ Name)) %>%
  select(seqid,start,end,Name) %>% 
  
  # Filter out genes that may have >0.25 reads or contigs alignment coverage

  left_join(
    EV_MAB_reads_coverage_geneID_GO_terms_synteny 
  ) %>%
  as.data.frame() %>% 
  filter(reads_used == "allMapped") %>%
  filter(mapping_reads_source  != "self_reads") %>%
  filter(fraction_overlap  < 0.25) %>%
  filter(syn_cov  < 0.25) %>%
  distinct() %>%
  # select(seqid,start,end,Name,fraction_overlap,syn_cov) %>%
  select(-overlap_count,-reads_used,-query,-reads_genome,-overlap_length, -query_genome,-total_length_ref,-ref_genome,-V4,-V5,-V6,-V7,-GO_terms) %>%
  distinct(seqid,start, end, genome, gene_id, fraction_overlap) 


#summary 
EV_MAB_specific_genes_nogo_noblastp %>%
  distinct(gene_id,genome) %>%
  group_by(genome) %>%
  summarise(count = n())




# read EVMAB genes that lack GO terms annotation

phylotools::read.fasta(paste0(DATA_DIR,"EV_MAB.nogo_noreads.prot.fasta")) %>%nrow()
  left_join(EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>%
              # distinct(qseqid) %>%
              rename(seq.name=qseqid)) %>%
  mutate(stitle_bp_annot = str_replace_na(stitle_bp_annot,"NA")) %>%
  filter(stitle_bp_annot == "NA") %>%
  select(-stitle_bp_annot,-genome) %>%
  phylotools::dat2fasta(paste0(DATA_DIR,"EV_MAB.prot.lack_goterms_blastp.annot.fasta"))


## Final EVMAB genes that lack GO-terms annotated but showed blastp hits agains NCBI-NR, uniprot-sprot and AT and OS proteoms 
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits.frac_overlap <- 
EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>% 
  separate(qseqid,into = c("Name","isomers"), sep = "-") %>%
  # filter(genome == "MA") %>%
  select(-genome) %>% 
  left_join(
    EV_MAB_reads_coverage_geneID %>%
      mutate(Name = case_when(str_detect(Name,"Ma") |
                                str_detect(Name,"Mb")~ str_c(Name,".1"),
                                      TRUE ~ Name))
    ) %>%
  as.data.frame() %>% 
  mutate(isomers = str_replace_na(isomers, "")) %>%
  mutate(Name = case_when(str_detect(Name,"-") ~ str_c(Name,isomers,sep = "-"),
                          TRUE ~ Name)) %>%
  filter(reads_used   == "allMapped") %>%
  mutate(
    mapping_reads_source = case_when(
      ref_genome == "ev_mazia" &
        reads_genome == "EV genotypes" ~ "self_reads",
      ref_genome == "ev_bedadeti" &
        reads_genome == "EV genotypes" ~ "self_reads",
      ref_genome == "musa_acuminata" &
        reads_genome == "A-sub genome banana genotypes" ~ "self_reads",
      ref_genome == "musa_balbisiana" &
        reads_genome == "B-sub genome banana genotypes" ~ "self_reads",
      TRUE ~ "cross_species_reads"
    )
  ) %>%
  filter(reads_used == "allMapped") %>%
  filter(mapping_reads_source  != "self_reads") %>%
  filter(fraction_overlap  < 0.25) %>%
  select(-overlap_count,-reads_used,-query,-reads_genome,-overlap_length,-total_length_ref,-group,-ref_genome) %>%
  mutate(fraction_overlap = round(fraction_overlap,2)) %>%
  distinct() 


## Final EVMAB genes that lack GO-terms annotated but showed blastp hits agains NCBI-NR, uniprot-sprot and AT and OS proteoms 

EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits.frac_overlap_syn <-
  EV_MAB_nogo_noreads_ncbi_uniprot_atos_hits %>% 
    separate(qseqid,into = c("Name","isomers"), sep = "-") %>%
    # filter(genome == "MA") %>%
    select(-genome) %>%
    left_join(
      EV_MAB_reads_coverage_geneID %>%
        mutate(Name = case_when(str_detect(Name,"Ma") |
                                  str_detect(Name,"Mb")~ str_c(Name,".1"),
                                TRUE ~ Name))
    ) %>%
    as.data.frame() %>% 
    mutate(isomers = str_replace_na(isomers, "")) %>%
    mutate(Name = case_when(str_detect(Name,"-") ~ str_c(Name,isomers,sep = "-"),
                            TRUE ~ Name)) %>%
    select(seqid,start,end,Name) %>%
    # head() %>%
    left_join(
      EV_MAB_reads_coverage_geneID_GO_terms_synteny 
    ) %>%
    as.data.frame() %>% 
    filter(reads_used == "allMapped") %>%
    filter(mapping_reads_source  != "self_reads") %>%
    filter(fraction_overlap  < 0.25) %>%
    filter(syn_cov  < 0.25) %>%
  distinct() %>%
    # select(seqid,start,end,Name,fraction_overlap,syn_cov) %>%
    select(-overlap_count,-reads_used,-query,-reads_genome,-overlap_length, -query_genome,-total_length_ref,-ref_genome,-V4,-V5,-V6,-V7,-gene_id,-GO_terms) %>%
  distinct()

   


# output their coding regions in bed format



