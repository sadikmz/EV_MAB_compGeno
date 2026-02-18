suppressPackageStartupMessages({
  library(tidyverse)
  library(phylotools)
})

# Exclude TE-inserted gene models from MAKER-predicted proteins.
# EV landraces (Mazia, Bedadeti) are filtered; Musa/E.glaucum retained as-is
# for comparative downstream analyses.

msg <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

TESORTER_DIR <- "tesorter_out"
FASTA_DIR    <- "gene_prediction_EDTA_masked"
EV_SPECIES   <- c("EV_mazia", "EV_bedadeti")
SPECIES_LIST <- c("EV_bedadeti", "EV_mazia", "Ensete_glaucum",
                  "Musa_acuminata", "Musa_balbisiana")

# Step 1: Load TEsorter domain (.dom.tsv) and classification (.cls.tsv) data
msg("Loading TEsorter results from: ", TESORTER_DIR)
load_te <- function(sp, ext) {
  path <- file.path(TESORTER_DIR, paste0(sp, "_TE.rexdb-plant.", ext, ".tsv"))
  tryCatch(read.delim(path) %>% mutate(genome = sp),
           error = function(e) stop(sprintf("Failed to read %s for %s: %s",
                                            ext, sp, conditionMessage(e))))
}
TEsorter_annotation     <- lapply(SPECIES_LIST, load_te, "dom") |> bind_rows()
TEsorter_annotation_cls <- lapply(SPECIES_LIST, load_te, "cls") |> bind_rows()
msg("Loaded: ", nrow(TEsorter_annotation), " dom rows; ",
    nrow(TEsorter_annotation_cls), " cls rows")

# Step 2: Reformat TEsorter annotations (MUST be defined before any join)
msg("Reformatting TEsorter annotations")
TEsorter_annotation.v1 <- TEsorter_annotation %>%
  mutate(
    seq.name = str_remove(X.id, "Class_I+.*[^|]"),
    seq.name = str_remove(seq.name, "\\|"),
    seq.name = str_remove(seq.name, ":\\d+-\\d+")
  ) %>%
  select(length, evalue, coverage, genome, seq.name) %>%
  left_join(
    TEsorter_annotation_cls %>% mutate(seq.name = str_remove(X.TE, ":\\d+-\\d+")),
    by = c("seq.name", "genome")
  ) %>%
  select(seq.name, genome, length, coverage, evalue, Order, Superfamily, Clade) %>%
  mutate(TE = str_c(Order, Superfamily, Clade, sep = ":"))
msg("TEsorter_annotation.v1 rows: ", nrow(TEsorter_annotation.v1))

# Step 3: Load predicted protein FASTAs for all species
msg("Loading protein FASTAs from: ", FASTA_DIR)
Ensete_MUSA_seq <- lapply(SPECIES_LIST, function(sp) {
  tryCatch(read.fasta(file.path(FASTA_DIR, paste0(sp, ".fasta"))),
           error = function(e) stop(sprintf("Failed to read FASTA for %s: %s",
                                            sp, conditionMessage(e)))) %>%
    mutate(
      genome   = sp,
      seq.name = str_extract(seq.name, "^\\S+")   # strip description text after first space
    )
}) |> bind_rows()
msg("Total sequences loaded: ", nrow(Ensete_MUSA_seq))

# Step 4: Exclude TE-inserted proteins for EV landraces; write filtered FASTAs
msg("Excluding TE-inserted gene models for EV landraces")
walk(EV_SPECIES, function(sp) {
  out <- file.path(FASTA_DIR, paste0(sp, ".TE_excluded.fasta"))
  te_excl <- Ensete_MUSA_seq %>%
    filter(genome == sp) %>%
    left_join(TEsorter_annotation.v1, by = c("seq.name", "genome")) %>%
    filter(is.na(length)) %>%
    select(seq.name, seq.text)
  dat2fasta(te_excl, out)
  msg("Written: ", out, " | TE-excluded proteins: ", nrow(te_excl))
})

# Step 5: Per-species summary table
msg("Generating TE exclusion summary")
summary_tbl <- Ensete_MUSA_seq %>%
  left_join(TEsorter_annotation.v1, by = c("seq.name", "genome")) %>%
  distinct(seq.name, genome, .keep_all = TRUE) %>%   # dedup multi-domain hits per protein
  mutate(has_TE = !is.na(Order)) %>%
  group_by(genome) %>%
  summarise(total = n(), te_inserted = sum(has_TE),
            te_excluded = total - te_inserted, .groups = "drop")
message(paste(capture.output(print(summary_tbl)), collapse = "\n"))

summary_out <- file.path(FASTA_DIR, "TE_exclusion_summary.tsv")
write_tsv(summary_tbl, summary_out)
msg("Summary written: ", summary_out)
