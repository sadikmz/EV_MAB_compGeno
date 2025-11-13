

interproscan_GO_terms.v1.MAB <- 
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "MA_specific") %>%
  mutate(
    description = str_replace(
      description,
      "isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",
      "isopentenyl diphosphate biosynthetic process..."
    )
  ) %>%
  dplyr::filter(description != "NA") %>%
  # dplyr::filter(Count >  1 ) %>%
  filter(!str_detect(description, "respiratory electron transport chain"),
         !str_detect(description, "phosphoenolpyruvate"),
         !str_detect(description, "photosynthesis"),
         !str_detect(description, "ATP metabolic process"),
         !str_detect(description, "cellular oxidant detoxification"),
         !str_detect(description, "microtubule-based movement"),
         !str_detect(description, "photosynthetic electron transport chain"),
         !str_detect(description, "vesicle-mediated transport"),
         !str_detect(description, "actin filament organization"),
         !str_detect(description, "nucleosome assembly"),
         !str_detect(description, "regulation of proton transport"),
         
         !str_detect(description, "translation"),
         # MA addtional
         !str_detect(description, "vacuolar transport"),
         !str_detect(description, "folding"),
         !str_detect(description, "prolyl-tRNA aminoacylation"),
         
         !str_detect(description, "RNA processing"),
         !str_detect(description, "protein ubiquitination"),
         !str_detect(description, "negative regulation of organ growth"),
         !str_detect(description, "mitochondrial")
  ) %>%
  select(category, description, Count) %>%
  mutate(category = str_remove(category,"_specific")) %>%
  rename(#cat = category,
         MA = Count) %>%
  full_join(

interproscan_GO_terms.v1 %>%
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "MB_specific") %>%
  dplyr::filter(description != "NA") %>%
  filter(!str_detect(description, "respiratory electron transport chain"),
         !str_detect(description, "phosphoenolpyruvate"),
         !str_detect(description, "photosynthesis"),
         !str_detect(description, "ATP metabolic process"),
         !str_detect(description, "cellular oxidant detoxification"),
         !str_detect(description, "microtubule-based movement"),
         !str_detect(description, "photosynthetic electron transport chain"),
         !str_detect(description, "vesicle-mediated transport"),
         !str_detect(description, "actin filament organization"),
         !str_detect(description, "nucleosome assembly"),
         !str_detect(description, "regulation of proton transport"),
         
         !str_detect(description, "translation"),
         # MA addtional
         !str_detect(description, "vacuolar transport"),
         !str_detect(description, "RNA processing"),
         !str_detect(description, "protein ubiquitination"),
         !str_detect(description, "negative regulation of organ growth"),
         !str_detect(description, "mitochondrial")) %>%
  select(category, description, Count) %>%
  mutate(category = str_remove(category,"_specific")) %>%
  rename(cat = category,
         MB = Count)) %>%
  select(description,MA,MB) %>%
  # head()
  mutate(MA = str_replace_na(MA, "0"),
         MB = str_replace_na(MB, "0"),
         MA = as.numeric(MA),
         MB = as.numeric(MB)) %>%
  # pivot_longer(cols = c(MA,MB),values_to = "count") %>%
  as.data.frame()



interproscan_GO_terms.v1.MAB %>%
  column_to_rownames("description")


## EV mazia and bedadeti 

interproscan_GO_terms.v1.EV <- 
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "EV_bedadeti_specific") %>%
  mutate(
    description = str_replace(description, "isopentenyl diphosphate biosynthetic process...", "isopentenyl diphosphate biosynthetic process"),
    # description = str_replace(description, "transmembrane transport", "transmembrane stransport")
  ) %>%
  # head()
  filter(!str_detect(description, "proteolysis"),
         description != "lysine biosynthetic process via diaminopimelate", 
         description != "iron-sulfur cluster assembly",
         description != "iron-sulfur cluster assembly",
         description != "transsulfuration",
         description != "amino acid metabolic process",
         description != "L-histidine biosynthetic process") %>%
  dplyr::filter(description != "NA") %>%
  dplyr::filter(Count >=  2 ) %>%
  # distinct(description)
  
  rbind(
    interproscan_GO_terms.v1 %>% 
      dplyr::rename(Count = total_genes,
                    category = clust) %>%
      dplyr::filter(category == "EV_bedadeti_specific")  ) %>%
  select(-category,-gene_list) %>%
  rename(EV_bedadeti = Count) %>%
  full_join(
interproscan_GO_terms.v1 %>% 
  dplyr::rename(Count = total_genes,
                category = clust) %>%
  dplyr::filter(category == "EV_mazia_specific") %>%
  select(-category,-gene_list) %>%
  rename(EV_mazia = Count) ) %>%
  mutate(
    description = str_replace(description, "intracellular signal transduction","Intracellular signal transduction"),
    description = str_replace(description, "regulation of DNA-templated transcription","Regulation of DNA-templated transcription"),
    description = str_replace(description, "signal transduction system","Signal transduction system"),
    description = str_replace(description, "carbohydrate metabolic process","Carbohydrate metabolic process"),
    description = str_replace(description, "proteasome assembly","Proteasome assembly"),
    description = str_replace(description, "glycolytic process","Glycolytic process"),
    description = str_replace(description, "proton transmembrane transport","Proton transmembrane transport"),
    description = str_replace(description, "lipid metabolic process","Lipid metabolic process"),
    description = str_replace(description, "biosynthetic process","Biosynthetic process"),
    description = str_replace(description, "protein phosphorylation","Protein phosphorylation"),
    description = str_replace(description, "transmembrane transport","Transmembrane transport"),
    description = str_replace(description, "regulation of DNA-templated transcription","Regulation of DNA-templated transcription"),
    description = str_replace(description, "phosphorelay signal transduction system","Phosphorelay signal transduction system"),
    description = str_replace(description, "phosphorylation","Phosphorylation"),
    description = str_replace(description, "signal transduction","Signal transduction"),
    description = str_replace(description, "carbohydrate metabolic process","Carbohydrate metabolic process"),
    description = str_replace(description, "amino acid metabolic process","Amino acid metabolic process"),
    description = str_replace(description, "glycolytic process","Glycolytic process"),
    description = str_replace(description, "gluconeogenesis","Gluconeogenesis"),
    description = str_replace(description, "glutathione metabolic process","Glutathione metabolic process"),
    description = str_replace(description, "transsulfuration","Transsulfuration"),
    description = str_replace(description, "protein folding","Protein folding")
  ) 

  