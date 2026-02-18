suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
  library(phylotools)
})

# --- Logging helper ---
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

# --- Parameters ---
TESORTER_DIR <- "tesorter_out"
FASTA_DIR    <- "gene_prediction_EDTA_masked"
PLOT_DIR     <- "tesorter_out"
SPECIES_LIST <- c("EV_bedadeti", "EV_mazia", "Ensete_glaucum",
                  "Musa_acuminata", "Musa_balbisiana")

GENOME_LABELS <- c(
  "Musa_balbisiana" = "MB\nproteins",
  "Musa_acuminata"  = "MA\nproteins",
  "EV_mazia"        = "EV (Mazia)\nproteins (AED<=0.25)",
  "Ensete_glaucum"  = "EG",
  "EV_bedadeti"     = "EV (Bedadeti)\nproteins (AED<=0.25)"
)

# Ordered TE family levels (Order:Superfamily:Clade)
TE_LEVELS <- c(
  "DIRS:unknown", "TIR:Tc1_Mariner", "TIR:Sola1", "TIR:PiggyBac",
  "TIR:PIF_Harbinger", "TIR:MuDR_Mutator", "TIR:hAT",
  "pararetrovirus:unknown", "mixture", "LINE:unknown", "LTR:mixture",
  "LTR:Copia:Bianca", "LTR:Copia", "LTR:Copia:Bryco", "LTR:Copia:Alesia",
  "LTR:Copia:Gymco-I", "LTR:Copia:Gymco-III", "LTR:Copia:Gymco-II",
  "LTR:Copia:Gymco-IV", "LTR:Copia:Lyco", "LTR:Copia:Osser",
  "LTR:Copia:TAR", "LTR:Copia:mixture", "LTR:Copia:Tork", "LTR:Copia:Ale",
  "LTR:Copia:Angela", "LTR:Copia:SIRE", "LTR:Copia:Ikeros",
  "LTR:Copia:Ivana", "LTR:Gypsy:Phygy", "LTR:Gypsy:Selgy",
  "LTR:Gypsy:Athila", "LTR:Gypsy:TatI", "LTR:Gypsy:chromo-outgroup",
  "LTR:Gypsy:chromo-unclass", "LTR:Gypsy:Chlamyvir", "LTR:Gypsy:TatII",
  "LTR:Gypsy:non-chromo-outgroup", "LTR:Gypsy:TatIII", "LTR:Gypsy:Tcn1",
  "LTR:Gypsy:Ogre", "LTR:Gypsy:mixture", "LTR:Gypsy:Tekay",
  "LTR:Gypsy:CRM", "LTR:Gypsy:Galadriel", "LTR:Gypsy:Reina",
  "LTR:Gypsy:Retand"
)

# Ordered TE family levels for the family-count plot (Superfamily:Clade only)
TE_FAM_LEVELS <- c(
  "DIRS", "Tc1_Mariner", "Sola1", "PiggyBac", "PIF_Harbinger",
  "MuDR_Mutator", "hAT", "pararetrovirus", "LINE", "mixture",
  "Bianca", "Copia", "Bryco", "Alesia", "Gymco-I", "Gymco-III",
  "Gymco-II", "Gymco-IV", "Lyco", "Osser", "TAR", "Tork", "Ale",
  "Angela", "SIRE", "Ikeros", "Ivana", "Phygy", "Selgy", "Athila",
  "TatI", "chromo-outgroup", "chromo-unclass", "Chlamyvir", "TatII",
  "non-chromo-outgroup", "TatIII", "Tcn1", "Ogre", "Tekay",
  "CRM", "Galadriel", "Reina", "Retand"
)

PLOT_COLOURS <- c("black", "blue", "green", "red", "gray", "orange", "#808000")

# Shared ggplot theme (bold axes/strips, no legend title)
base_theme <- theme_bw() + theme(
  axis.title.x          = element_text(face = "bold"),
  axis.title.y          = element_text(face = "bold"),
  axis.text.x           = element_text(face = "bold"),
  axis.text.y           = element_text(face = "bold"),
  strip.text.x          = element_text(face = "bold"),
  legend.text           = element_markdown(margin = margin(r = 10), face = "bold"),
  legend.title          = element_blank(),
  legend.box.background = element_rect(colour = "black")
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Save a ggplot to both TIFF and PNG
save_plot <- function(p, out_base, width, height) {
  tiff_path <- paste0(out_base, ".tiff")
  png_path  <- paste0(out_base, ".png")
  ggsave(tiff_path, plot = p, width = width, height = height,
         device = "tiff", dpi = 300, compression = "lzw")
  msg("Saved: ", tiff_path)
  ggsave(png_path,  plot = p, width = width, height = height,
         device = png, type = "cairo", dpi = 300)
  msg("Saved: ", png_path)
  invisible(p)
}

#' Load TEsorter dom.tsv and cls.tsv for all species
load_tesorter_data <- function(species_list, base_dir) {
  read_sp <- function(sp, suffix) {
    path <- file.path(base_dir, paste0(sp, "_TE.rexdb-plant.", suffix, ".tsv"))
    tryCatch(
      read.delim(path) %>% mutate(genome = sp),
      error = function(e) stop(sprintf(
        "Cannot read %s for %s (%s): %s", suffix, sp, path, conditionMessage(e)))
    )
  }
  dom <- map_dfr(species_list, read_sp, suffix = "dom")
  cls <- map_dfr(species_list, read_sp, suffix = "cls")
  msg("  Loaded dom rows: ", nrow(dom), " | cls rows: ", nrow(cls))
  list(annotation = dom, annotation_cls = cls)
}

# ---------------------------------------------------------------------------
# Step 1: Load TEsorter data
# ---------------------------------------------------------------------------
msg("=== parse_plot_TEsorter_out ===")
msg("Loading TEsorter data from: ", TESORTER_DIR)

tesorter_data           <- load_tesorter_data(SPECIES_LIST, TESORTER_DIR)
TEsorter_annotation     <- tesorter_data$annotation
TEsorter_annotation_cls <- tesorter_data$annotation_cls

# ---------------------------------------------------------------------------
# Step 2: Reformat TEsorter annotations
# ---------------------------------------------------------------------------
msg("Reformatting TEsorter annotations into v1 table")

TEsorter_annotation.v1 <-
  TEsorter_annotation %>%
  mutate(
    seq.name = str_remove(X.id, "Class_I+.*[^|]"),
    seq.name = str_remove(seq.name, "\\|"),
    seq.name = str_remove(seq.name, ":\\d+-\\d+")
  ) %>%
  select(length, evalue, coverage, genome, seq.name) %>%
  left_join(
    TEsorter_annotation_cls %>%
      mutate(seq.name = str_remove(X.TE, ":\\d+-\\d+")),
    by = c("seq.name", "genome")
  ) %>%
  select(seq.name, genome, length, coverage, evalue, Order, Superfamily, Clade) %>%
  mutate(TE = str_c(Order, Superfamily, Clade, sep = ":"))

msg("  TEsorter_annotation.v1 rows: ", nrow(TEsorter_annotation.v1))

# ---------------------------------------------------------------------------
# Step 3: Load predicted protein FASTAs
# ---------------------------------------------------------------------------
msg("Loading predicted protein FASTAs from: ", FASTA_DIR)

Ensete_MUSA_seq <- map_dfr(SPECIES_LIST, function(sp) {
  path <- file.path(FASTA_DIR, paste0(sp, ".fasta"))
  tryCatch(
    read.fasta(path),
    error = function(e) stop(sprintf(
      "Cannot read FASTA for %s (%s): %s", sp, path, conditionMessage(e)))
  ) %>%
    mutate(
      genome   = sp,
      seq.name = str_extract(seq.name, "^\\S+")   # strip description text after first space
    )
})

msg("Total sequences loaded: ", nrow(Ensete_MUSA_seq))

# ---------------------------------------------------------------------------
# Step 4: Join predicted proteins with TEsorter classifications
# ---------------------------------------------------------------------------
msg("Joining predicted proteins with TEsorter classifications")

predicted_prot_TEsorter_joined.AED <-
  Ensete_MUSA_seq %>%
  mutate(seq.len = nchar(seq.text)) %>%
  select(-seq.text) %>%
  left_join(TEsorter_annotation.v1, by = c("seq.name", "genome")) %>%
  mutate(
    coverage = round(length / seq.len, 2),
    TE       = str_replace_na(TE, "No_TE_found")
  ) %>%
  filter(TE != "No_TE_found")

msg("  Joined TE-annotated proteins: ", nrow(predicted_prot_TEsorter_joined.AED))

# ---------------------------------------------------------------------------
# Step 5: Build v2 table with cleaned labels and factor levels
# ---------------------------------------------------------------------------
TEsorter_annotation.v2 <-
  predicted_prot_TEsorter_joined.AED %>%
  mutate(
    TE       = str_replace(TE, "mixture\\:mixture", "mixture"),
    TE       = str_remove(TE, "\\:unknown"),
    genome   = str_replace_all(genome, GENOME_LABELS),
    coverage = case_when(coverage > 1 ~ 1, TRUE ~ coverage),
    genome   = factor(genome, levels = c(
      "EV (Mazia)\nproteins (AED<=0.25)",
      "EV (Bedadeti)\nproteins (AED<=0.25)",
      "EG", "MA\nproteins", "MB\nproteins"
    )),
    TE = factor(TE, levels = TE_LEVELS)
  )

# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

#' TE domain coverage across predicted gene models
plot_te_coverage <- function(data, out_base) {
  # "EG" = Ensete glaucum; "MS" retained from original filter — verify if "MB" was intended
  p <- data %>%
    filter(!genome %in% c("EG", "MS")) %>%
    ggplot(aes(TE, coverage, colour = genome)) +
    geom_jitter(width = 0.1, height = 0.1, size = 3, alpha = 0.5) +
    facet_grid(~genome) +
    coord_flip() +
    scale_colour_manual(values = PLOT_COLOURS) +
    labs(
      y = "Length coverage (%) of TE-inserted genes in the predicted proteins",
      x = "Transposable elements\nOrder:Superfamily:Clade"
    ) +
    base_theme +
    theme(
      axis.title.x = element_text(size = 9,  vjust = 1),
      axis.title.y = element_text(size = 9,  vjust = -8),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      strip.text.x = element_text(size = 9),
      legend.text  = element_markdown(margin = margin(r = 10), size = 10),
      legend.position = "none"
    )
  save_plot(p, out_base, width = 8, height = 6)
}

#' Number of predicted proteins with TE insertions per TE family
plot_te_family_counts <- function(data, out_base) {
  # "EG" = Ensete glaucum; "MS" retained from original filter — verify if "MB" was intended
  # str_detect(genome, "^M") additionally excludes Musa species whose labels start with "M"
  p <- data %>%
    filter(!genome %in% c("EG", "MS"), !str_detect(genome, "^M"),
           Order != "mixture") %>%
    group_by(genome, TE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      TE     = str_remove(str_remove(TE, "LTR:"), ":unknown"),
      type   = str_extract(TE, "\\w+:"),
      TE_fam = str_remove(TE, "\\w+:"),
      type   = str_replace_na(type, "NA")
    ) %>%
    filter(TE != "mixture") %>%
    { bind_rows(
        filter(., type == "NA") %>% mutate(type = TE),
        filter(., type != "NA")
      ) } %>%
    mutate(
      type   = str_remove(type, ":"),
      TE_fam = factor(TE_fam, levels = TE_FAM_LEVELS)
    ) %>%
    ggplot(aes(TE_fam, count, colour = type)) +
    geom_point(size = 4, alpha = 0.7) +
    facet_grid(~genome) +
    coord_flip() +
    # Break values are data-derived: observed quantiles of TE-family counts per genome
    scale_y_continuous(breaks = c(50, 2500, 5000, 6746, 8651)) +
    labs(
      y = "Number of TE-encoding protein families inserted in the predicted proteins",
      x = "Transposable elements\nSuperfamily:Clade"
    ) +
    base_theme +
    theme(
      axis.title.x    = element_text(size = 9,  vjust = 1),
      axis.title.y    = element_text(size = 9,  vjust = -9),
      axis.text.x     = element_text(size = 8),
      axis.text.y     = element_text(size = 9),
      strip.text.x    = element_text(size = 9),
      legend.text     = element_markdown(margin = margin(r = 10), size = 10),
      legend.position = c(0.85, 0.8),
      legend.key.size = unit(0.1, "cm")
    )
  save_plot(p, out_base, width = 7, height = 6)
}

#' Frequency of distinct TE classes per predicted gene model
plot_te_frequency <- function(data, out_base) {
  p <- data %>%
    group_by(genome, seq.name) %>%
    summarise(TE_count = n(), .groups = "drop") %>%
    group_by(genome, TE_classes_per_gene = TE_count) %>%
    summarise(n_gene_models = n(), .groups = "drop") %>%
    filter(n_gene_models != 0) %>%
    ggplot(aes(TE_classes_per_gene, n_gene_models, colour = genome)) +
    geom_point(size = 4) +
    scale_colour_manual(values = PLOT_COLOURS) +
    # Break values are data-derived: observed gene-model count landmarks across genomes
    scale_y_continuous(breaks = c(1, 500, 2500, 4500, 6500, 9000)) +
    labs(
      x = "Number of distinct TE classes inserted per single protein-coding gene",
      y = "Number of TE-containing gene models"
    ) +
    base_theme +
    theme(
      axis.title.x    = element_text(size = 16),
      axis.title.y    = element_text(size = 16),
      axis.text.x     = element_text(size = 12),
      axis.text.y     = element_text(size = 12),
      strip.text.x    = element_text(size = 12),
      legend.text     = element_markdown(margin = margin(r = 10), size = 12),
      legend.position = "top"
    )
  save_plot(p, out_base, width = 5, height = 5)
}

# ---------------------------------------------------------------------------
# Step 6: Generate all figures
# ---------------------------------------------------------------------------
msg("Generating figures")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_te_coverage(
  TEsorter_annotation.v2,
  file.path(PLOT_DIR, "TE_sorter.EV_AED25.MAB")
)

plot_te_family_counts(
  TEsorter_annotation.v2,
  file.path(PLOT_DIR, "TE_sorter.EV_AED25.MAB.TE_families_count")
)

plot_te_frequency(
  TEsorter_annotation.v2,
  file.path(PLOT_DIR, "TE_sorter.EV_MAB.frequency")
)

# ---------------------------------------------------------------------------
# Step 7: Summary table
# ---------------------------------------------------------------------------
msg("Writing summary table")

summary_out <- file.path(PLOT_DIR, "count_TE_inserted_gene.txt")

TEsorter_annotation.v2 %>%
  group_by(genome, seq.name) %>%
  summarise(TE_count = n(), .groups = "drop") %>%
  group_by(genome) %>%
  summarise(count = n(), max = max(TE_count), .groups = "drop") %>%
  write.table(summary_out, col.names = TRUE, row.names = FALSE,
              quote = FALSE, sep = "\t")

msg("Summary table written: ", summary_out)
msg("=== parse_plot_TEsorter_out complete ===")
