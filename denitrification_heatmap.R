#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# USER SETTINGS
# =========================================================

orftable_file <- "/hpcfs/home/w20012148/SqueezeMeta/1st_run/results/13.1st_run.orftable"
kegg_diamond_file <- "/hpcfs/home/w20012148/SqueezeMeta/1st_run/intermediate/04.1st_run.kegg.diamond"

out_dir <- "/hpcfs/home/w20012148/SqueezeMeta/1st_run/results/paper_results/fig_3/now_with_closness/results"

samples_to_use <- c("Sample21", "Sample22", "Sample23")

use_log10_for_colour <- TRUE
plot_pseudocount     <- 1e-6

png_width  <- 7000
png_height <- 5200
png_res    <- 450

# =========================================================
# CONFIDENCE THRESHOLDS
# Same thresholds as nitrification / assimilatory figures
# =========================================================

high_identity_threshold       <- 70
high_query_coverage_threshold <- 70

medium_identity_threshold       <- 40
medium_query_coverage_threshold <- 40

acceptable_evalue_threshold <- 1e-5
minimum_alignment_length_aa <- 80

# =========================================================
# DENITRIFICATION DEFINITIONS
# Multiple acceptable KEGG IDs are separated by semicolons.
# =========================================================

denit_defs <- tibble::tribble(
  ~gene,  ~display, ~kegg_ids,                       ~gene_regex,
  "narG", "narG",   "K00370",                       "\\bnarG\\b|nitrate reductase alpha",
  "narH", "narH",   "K00371",                       "\\bnarH\\b|nitrate reductase beta",
  "narI", "narI",   "K00374",                       "\\bnarI\\b|nitrate reductase gamma",
  "napA", "napA",   "K02567",                       "\\bnapA\\b|periplasmic nitrate reductase",
  "napB", "napB",   "K02568",                       "\\bnapB\\b|periplasmic nitrate reductase.*beta",
  "nirK", "nirK",   "K00368",                       "\\bnirK\\b|nitrite reductase.*copper",
  "nirS", "nirS",   "K15864",                       "\\bnirS\\b|cytochrome cd1 nitrite reductase",
  "norB", "norB",   "K04561", "\\bnorB\\b|nitric oxide reductase.*subunit B",
  "norC", "norC",   "K02305", "\\bnorC\\b|nitric oxide reductase.*subunit C",
  "nosZ", "nosZ",   "K00376",                       "\\bnosZ\\b|nitrous-oxide reductase"
)

gene_order <- denit_defs$display

gene_to_reaction <- c(
  "narG" = "Nitrate → Nitrite",
  "narH" = "Nitrate → Nitrite",
  "narI" = "Nitrate → Nitrite",
  "napA" = "Nitrate → Nitrite",
  "napB" = "Nitrate → Nitrite",
  "nirK" = "Nitrite → NO",
  "nirS" = "Nitrite → NO",
  "norB" = "NO → N2O",
  "norC" = "NO → N2O",
  "nosZ" = "N2O → N2"
)

reaction_levels <- c(
  "Nitrate → Nitrite",
  "Nitrite → NO",
  "NO → N2O",
  "N2O → N2"
)

# =========================================================
# HELPERS
# =========================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

extract_rank <- function(tax_string, prefix) {
  out <- str_match(tax_string, paste0("(^|;)", prefix, "([^;]+)"))[, 3]
  out[is.na(out)] <- ""
  out
}

clean_tax_label <- function(x) {
  x <- str_replace_all(x, "^Candidatus\\s+", "Ca. ")
  x <- str_replace_all(x, "_", " ")
  trimws(x)
}

is_bad_exact <- function(x) {
  tolower(trimws(x)) %in% c(
    "", "unclassified", "unclassified bacteria", "bacteria",
    "unknown", "uncultured bacterium"
  )
}

normalise_colname <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

first_existing_col <- function(df, candidates) {
  nms <- colnames(df)
  norm_nms <- normalise_colname(nms)
  norm_candidates <- normalise_colname(candidates)

  hit <- which(norm_nms %in% norm_candidates)
  if (length(hit) == 0) return(NA_character_)

  nms[hit[1]]
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", as.character(x))))
}

weighted_mean_safe <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)

  keep <- !is.na(x) & !is.na(w)
  if (!any(keep)) return(NA_real_)

  x <- x[keep]
  w <- w[keep]

  if (sum(w, na.rm = TRUE) > 0) {
    return(weighted.mean(x, w, na.rm = TRUE))
  }

  mean(x, na.rm = TRUE)
}

split_kegg_ids <- function(ids) {
  ids_vec <- trimws(unlist(strsplit(ids, ";", fixed = TRUE)))
  ids_vec[ids_vec != ""]
}

has_any_expected_ko <- function(kegg_id_col, expected_ids) {
  ids_vec <- split_kegg_ids(expected_ids)

  if (length(ids_vec) == 0) {
    return(rep(FALSE, length(kegg_id_col)))
  }

  Reduce(`|`, lapply(ids_vec, function(id) {
    str_detect(kegg_id_col, fixed(id))
  }))
}

match_kegg_gene <- function(kegg_id_col, gene_name_col, keggfun_col,
                            expected_kegg_ids, gene_regex) {

  has_expected_ko <- has_any_expected_ko(
    kegg_id_col,
    expected_kegg_ids
  )

  has_text_match <- str_detect(
    tolower(paste(gene_name_col, keggfun_col)),
    regex(tolower(gene_regex), ignore_case = TRUE)
  )

  has_expected_ko | has_text_match
}

match_kegg_hits_for_df <- function(df_in, defs_tbl) {
  long_list <- list()

  for (i in seq_len(nrow(defs_tbl))) {
    disp <- defs_tbl$display[i]
    expected_ids <- defs_tbl$kegg_ids[i]

    kegg_hits <- match_kegg_gene(
      kegg_id_col      = df_in$`KEGG ID`,
      gene_name_col    = df_in$`Gene name`,
      keggfun_col      = df_in$KEGGFUN,
      expected_kegg_ids = expected_ids,
      gene_regex       = defs_tbl$gene_regex[i]
    )

    if (any(kegg_hits)) {
      long_list[[length(long_list) + 1]] <- df_in[kegg_hits, ] %>%
        mutate(
          database = "KEGG",
          gene = disp,
          expected_kegg_ids = expected_ids
        )
    }
  }

  if (length(long_list) == 0) return(NULL)
  bind_rows(long_list)
}

read_kegg_diamond <- function(kegg_diamond_file) {

  if (!file.exists(kegg_diamond_file)) {
    stop("KEGG DIAMOND file not found: ", kegg_diamond_file)
  }

  message("Reading KEGG DIAMOND file: ", kegg_diamond_file)

  kegg_diam <- fread(
    kegg_diamond_file,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    data.table = FALSE
  )

  if (ncol(kegg_diam) != 12) {
    stop(
      "Expected 12 columns in KEGG DIAMOND file, but found ",
      ncol(kegg_diam),
      "."
    )
  }

  colnames(kegg_diam) <- c(
    "qseqid",
    "qlen",
    "sseqid",
    "slen",
    "pident",
    "alignment_length",
    "evalue",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send"
  )

  numeric_cols <- c(
    "qlen", "slen", "pident", "alignment_length",
    "evalue", "bitscore", "qstart", "qend", "sstart", "send"
  )

  for (cc in numeric_cols) {
    kegg_diam[[cc]] <- as_numeric_clean(kegg_diam[[cc]])
  }

  kegg_diam %>%
    mutate(
      orf_id_for_validation = as.character(qseqid),
      best_kegg_reference = as.character(sseqid),
      diamond_ko = str_extract(sseqid, "K\\d{5}"),

      query_aligned_span = abs(qend - qstart) + 1,
      target_aligned_span = abs(send - sstart) + 1,

      kegg_query_coverage = ifelse(
        !is.na(qlen) & qlen > 0,
        pmin((query_aligned_span / qlen) * 100, 100),
        NA_real_
      ),

      kegg_target_coverage = ifelse(
        !is.na(slen) & slen > 0,
        pmin((target_aligned_span / slen) * 100, 100),
        NA_real_
      )
    ) %>%
    transmute(
      orf_id_for_validation,
      best_kegg_reference,
      diamond_ko,
      kegg_pident = pident,
      kegg_alignment_length = alignment_length,
      kegg_query_coverage,
      kegg_target_coverage,
      kegg_evalue = evalue,
      kegg_bitscore = bitscore,
      qlen,
      slen,
      qstart,
      qend,
      sstart,
      send
    )
}

ko_is_expected <- function(diamond_ko, expected_kegg_ids) {
  if (is.na(diamond_ko) || is.na(expected_kegg_ids)) return(FALSE)
  diamond_ko %in% split_kegg_ids(expected_kegg_ids)
}

classify_kegg_confidence <- function(
    expected_kegg_ids,
    diamond_ko,
    pident,
    query_coverage,
    alignment_length,
    evalue
) {

  has_expected_ko <- ko_is_expected(
    diamond_ko,
    expected_kegg_ids
  )

  if (!has_expected_ko) {
    return("Unvalidated")
  }

  has_identity <- !is.na(pident)
  has_qcov     <- !is.na(query_coverage)
  has_evalue   <- !is.na(evalue)
  has_aln_len  <- !is.na(alignment_length)

  good_evalue <- !has_evalue ||
    evalue <= acceptable_evalue_threshold

  good_length <- !has_aln_len ||
    alignment_length >= minimum_alignment_length_aa

  if (
    has_identity &&
    has_qcov &&
    pident >= high_identity_threshold &&
    query_coverage >= high_query_coverage_threshold &&
    good_evalue &&
    good_length
  ) {
    return("High")
  }

  if (
    has_identity &&
    pident >= medium_identity_threshold &&
    (!has_qcov ||
       query_coverage >= medium_query_coverage_threshold) &&
    good_evalue &&
    good_length
  ) {
    return("Medium")
  }

  "Low"
}

confidence_to_score <- function(confidence) {
  dplyr::case_when(
    confidence == "High" ~ 3,
    confidence == "Medium" ~ 2,
    confidence == "Low" ~ 1,
    confidence == "Unvalidated" ~ 0,
    TRUE ~ NA_real_
  )
}

score_to_average_confidence <- function(score) {
  dplyr::case_when(
    is.na(score) ~ "No detected ORFs",
    score >= 2.5 ~ "High",
    score >= 1.5 ~ "Medium",
    score > 0 ~ "Low",
    TRUE ~ "Unvalidated"
  )
}

confidence_label <- function(category) {
  dplyr::case_when(
    category == "High" ~ "H",
    category == "Medium" ~ "M",
    category == "Low" ~ "L",
    category == "Unvalidated" ~ "U",
    TRUE ~ ""
  )
}

# =========================================================
# READ ORF TABLE
# =========================================================

message("Reading ORF table...")

dt <- fread(
  orftable_file,
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = "",
  comment.char = "#",
  check.names = FALSE,
  data.table = FALSE
)

colnames(dt) <- trimws(colnames(dt))

required_cols <- c("Tax", "KEGG ID", "KEGGFUN")
missing_cols <- setdiff(required_cols, colnames(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in ORF table: ",
    paste(missing_cols, collapse = ", ")
  )
}

tpm_cols <- paste0("TPM ", samples_to_use)
missing_tpm <- setdiff(tpm_cols, colnames(dt))

if (length(missing_tpm) > 0) {
  stop(
    "Missing TPM columns in ORF table: ",
    paste(missing_tpm, collapse = ", ")
  )
}

for (cc in c(
  "Tax",
  "Gene name",
  "KEGG ID",
  "KEGGFUN",
  "KEGGPATH"
)) {
  if (!cc %in% colnames(dt)) {
    dt[[cc]] <- ""
  }

  dt[[cc]][is.na(dt[[cc]])] <- ""
}

orf_id_col <- first_existing_col(
  dt,
  c(
    "ORF",
    "orf",
    "ORF ID",
    "orf_id",
    "query",
    "qseqid",
    "gene_id",
    "Gene ID"
  )
)

if (is.na(orf_id_col)) {
  warning(
    "No obvious ORF ID column found. Using first column as ORF ID: ",
    colnames(dt)[1]
  )

  orf_id_col <- colnames(dt)[1]
}

message(
  "Using ORF ID column for KEGG validation: ",
  orf_id_col
)

dt$orf_id_for_validation <- as.character(dt[[orf_id_col]])

dt$mean_tpm <- rowMeans(
  as.data.frame(
    lapply(
      dt[, tpm_cols, drop = FALSE],
      as.numeric
    )
  ),
  na.rm = TRUE
)

dt$mean_tpm[is.na(dt$mean_tpm)] <- 0

# =========================================================
# TAXONOMY PARSING
# =========================================================

dt$phylum <- clean_tax_label(
  extract_rank(dt$Tax, "p_")
)

dt$genus <- clean_tax_label(
  extract_rank(dt$Tax, "g_")
)

dt <- dt %>%
  mutate(
    phylum = ifelse(
      is_bad_exact(phylum),
      "",
      phylum
    ),
    genus = ifelse(
      is_bad_exact(genus),
      "",
      genus
    )
  ) %>%
  filter(phylum != "")

dt_genus <- dt %>%
  filter(genus != "") %>%
  mutate(
    group_level = "genus",
    group_name = genus,
    row_label = genus,
    row_type = "genus"
  )

dt_phylum <- dt %>%
  mutate(
    group_level = "phylum",
    group_name = phylum,
    row_label = phylum,
    row_type = "phylum_summary"
  )

# =========================================================
# MATCH KEGG DENITRIFICATION GENES
# =========================================================

message("Matching KEGG denitrification genes...")

hits_genus <- match_kegg_hits_for_df(
  dt_genus,
  denit_defs
)

hits_phylum <- match_kegg_hits_for_df(
  dt_phylum,
  denit_defs
)

if (
  is.null(hits_genus) &&
  is.null(hits_phylum)
) {
  stop(
    "No KEGG denitrification matches found."
  )
}

hits_long <- bind_rows(
  hits_genus,
  hits_phylum
)

phyla_kept <- hits_long %>%
  distinct(phylum) %>%
  pull(phylum)

dt_genus <- dt_genus %>%
  filter(phylum %in% phyla_kept)

dt_phylum <- dt_phylum %>%
  filter(phylum %in% phyla_kept)

hits_long <- hits_long %>%
  filter(phylum %in% phyla_kept)

# =========================================================
# READ KEGG DIAMOND AND VALIDATE MATCHES
# =========================================================

message(
  "Building KEGG denitrification validation table..."
)

kegg_diamond_all <- read_kegg_diamond(
  kegg_diamond_file
)

unique_hit_validation <- hits_long %>%
  mutate(
    row_type_priority = ifelse(
      row_type == "genus",
      1,
      2
    )
  ) %>%
  arrange(row_type_priority) %>%
  group_by(
    orf_id_for_validation,
    gene,
    expected_kegg_ids
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(-row_type_priority)

# Find all DIAMOND hits for candidate ORFs, then retain only hits
# matching one of the acceptable KOs for that denitrification gene.
candidate_diamond_hits <- unique_hit_validation %>%
  select(
    orf_id_for_validation,
    gene,
    expected_kegg_ids
  ) %>%
  inner_join(
    kegg_diamond_all,
    by = "orf_id_for_validation"
  ) %>%
  rowwise() %>%
  mutate(
    expected_ko_match = ko_is_expected(
      diamond_ko,
      expected_kegg_ids
    )
  ) %>%
  ungroup() %>%
  filter(expected_ko_match) %>%
  arrange(
    orf_id_for_validation,
    gene,
    desc(
      replace_na(
        kegg_bitscore,
        -Inf
      )
    ),
    replace_na(
      kegg_evalue,
      Inf
    )
  ) %>%
  group_by(
    orf_id_for_validation,
    gene,
    expected_kegg_ids
  ) %>%
  slice(1) %>%
  ungroup()

unique_hit_validation <- unique_hit_validation %>%
  left_join(
    candidate_diamond_hits,
    by = c(
      "orf_id_for_validation",
      "gene",
      "expected_kegg_ids"
    )
  )

unique_hit_validation$confidence <- mapply(
  classify_kegg_confidence,
  expected_kegg_ids =
    unique_hit_validation$expected_kegg_ids,
  diamond_ko =
    unique_hit_validation$diamond_ko,
  pident =
    unique_hit_validation$kegg_pident,
  query_coverage =
    unique_hit_validation$kegg_query_coverage,
  alignment_length =
    unique_hit_validation$kegg_alignment_length,
  evalue =
    unique_hit_validation$kegg_evalue
)

unique_hit_validation <- unique_hit_validation %>%
  rowwise() %>%
  mutate(
    confidence_score =
      confidence_to_score(confidence),

    exact_expected_ko_in_orftable =
      any(
        vapply(
          split_kegg_ids(
            expected_kegg_ids
          ),
          function(id) {
            str_detect(
              `KEGG ID`,
              fixed(id)
            )
          },
          logical(1)
        )
      )
  ) %>%
  ungroup()


# =========================================================
# RESOLVE COMPETING GENE ASSIGNMENTS
# One ORF can contribute to only one denitrification gene.
# Prefer the strongest KO-supported DIAMOND assignment.
# =========================================================

unique_hit_validation <- unique_hit_validation %>%
  arrange(
    orf_id_for_validation,
    desc(exact_expected_ko_in_orftable),
    desc(confidence_score),
    desc(replace_na(kegg_bitscore, -Inf)),
    desc(replace_na(kegg_pident, -Inf)),
    desc(replace_na(kegg_query_coverage, -Inf)),
    replace_na(kegg_evalue, Inf)
  ) %>%
  group_by(orf_id_for_validation) %>%
  slice(1) %>%
  ungroup()

# Remove losing alternative gene calls before abundance aggregation
# so an ORF cannot contribute TPM to more than one gene.
hits_long <- hits_long %>%
  semi_join(
    unique_hit_validation %>%
      distinct(orf_id_for_validation, gene, expected_kegg_ids),
    by = c("orf_id_for_validation", "gene", "expected_kegg_ids")
  )

# =========================================================
# ORF-LEVEL VALIDATION TABLE
# =========================================================

orf_validation_table <- unique_hit_validation %>%
  transmute(
    orf_id = orf_id_for_validation,
    gene,
    reaction = gene_to_reaction[gene],
    expected_kegg_ids,
    diamond_ko,
    exact_expected_ko_in_orftable,
    best_kegg_reference,
    kegg_pident,
    kegg_alignment_length,
    kegg_query_coverage,
    kegg_target_coverage,
    kegg_evalue,
    kegg_bitscore,
    confidence,
    confidence_score,
    phylum,
    genus,
    taxon = ifelse(
      genus != "",
      genus,
      phylum
    ),
    mean_tpm,
    KEGG_ID = `KEGG ID`,
    KEGGFUN,
    Tax
  ) %>%
  arrange(
    gene,
    factor(
      confidence,
      levels = c(
        "High",
        "Medium",
        "Low",
        "Unvalidated"
      )
    ),
    desc(
      replace_na(
        kegg_pident,
        -Inf
      )
    ),
    desc(
      replace_na(
        kegg_query_coverage,
        -Inf
      )
    ),
    replace_na(
      kegg_evalue,
      Inf
    )
  )

write.table(
  orf_validation_table,
  file.path(
    out_dir,
    "kegg_denitrification_orf_validation_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation_summary_by_gene <- orf_validation_table %>%
  count(
    gene,
    confidence,
    name = "n_orfs"
  ) %>%
  pivot_wider(
    names_from = confidence,
    values_from = n_orfs,
    values_fill = 0
  )

write.table(
  validation_summary_by_gene,
  file.path(
    out_dir,
    "kegg_denitrification_orf_validation_summary_by_gene.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# JOIN VALIDATION BACK TO TAXON-LEVEL HITS
# =========================================================

hits_long <- hits_long %>%
  left_join(
    unique_hit_validation %>%
      select(
        orf_id_for_validation,
        gene,
        expected_kegg_ids,
        diamond_ko,
        best_kegg_reference,
        kegg_pident,
        kegg_alignment_length,
        kegg_query_coverage,
        kegg_target_coverage,
        kegg_evalue,
        kegg_bitscore,
        confidence,
        confidence_score
      ),
    by = c(
      "orf_id_for_validation",
      "gene",
      "expected_kegg_ids"
    )
  ) %>%
  mutate(
    confidence =
      replace_na(
        confidence,
        "Unvalidated"
      ),
    confidence_score =
      replace_na(
        confidence_score,
        0
      )
  )

# =========================================================
# EXCLUDE UNVALIDATED KEGG ASSIGNMENTS FROM MAIN FIGURE
# Unvalidated calls remain in the ORF validation output but
# do not contribute to heatmap abundance or confidence.
# =========================================================

hits_long <- hits_long %>%
  filter(confidence != "Unvalidated")

# =========================================================
# TOTAL TPM FOR NORMALISATION
# =========================================================

genus_totals <- dt_genus %>%
  group_by(
    phylum,
    group_level,
    group_name,
    row_label,
    row_type
  ) %>%
  summarise(
    total_group_tpm =
      sum(
        mean_tpm,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

phylum_totals <- dt_phylum %>%
  group_by(
    phylum,
    group_level,
    group_name,
    row_label,
    row_type
  ) %>%
  summarise(
    total_group_tpm =
      sum(
        mean_tpm,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

group_totals <- bind_rows(
  genus_totals,
  phylum_totals
)

# =========================================================
# AGGREGATE KEGG ABUNDANCE AND CONFIDENCE
# =========================================================

agg <- hits_long %>%
  group_by(
    phylum,
    group_level,
    group_name,
    row_label,
    row_type,
    gene
  ) %>%
  summarise(
    gene_tpm =
      sum(
        mean_tpm,
        na.rm = TRUE
      ),

    n_orfs =
      n_distinct(
        orf_id_for_validation
      ),

    n_high =
      sum(
        confidence == "High",
        na.rm = TRUE
      ),

    n_medium =
      sum(
        confidence == "Medium",
        na.rm = TRUE
      ),

    n_low =
      sum(
        confidence == "Low",
        na.rm = TRUE
      ),

    n_unvalidated =
      sum(
        confidence == "Unvalidated",
        na.rm = TRUE
      ),

    weighted_confidence_score =
      weighted_mean_safe(
        confidence_score,
        mean_tpm
      ),

    mean_confidence_score =
      mean(
        confidence_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  left_join(
    group_totals,
    by = c(
      "phylum",
      "group_level",
      "group_name",
      "row_label",
      "row_type"
    )
  ) %>%
  mutate(
    norm_value = ifelse(
      total_group_tpm > 0,
      gene_tpm / total_group_tpm,
      0
    ),

    average_confidence =
      score_to_average_confidence(
        weighted_confidence_score
      ),

    confidence_label =
      confidence_label(
        average_confidence
      )
  )

rows_kept <- agg %>%
  distinct(
    phylum,
    group_level,
    group_name,
    row_label,
    row_type
  )

phylum_order <- agg %>%
  group_by(phylum) %>%
  summarise(
    phylum_signal =
      sum(
        gene_tpm,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(phylum_signal)
  ) %>%
  pull(phylum)

phylum_summary_order <- agg %>%
  filter(
    row_type ==
      "phylum_summary"
  ) %>%
  group_by(
    phylum,
    row_label,
    row_type
  ) %>%
  summarise(
    sort_signal =
      sum(
        norm_value,
        na.rm = TRUE
      ),
    total_gene_tpm =
      sum(
        gene_tpm,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

genus_order <- agg %>%
  filter(
    row_type ==
      "genus"
  ) %>%
  group_by(
    phylum,
    row_label,
    row_type
  ) %>%
  summarise(
    sort_signal =
      sum(
        norm_value,
        na.rm = TRUE
      ),
    total_gene_tpm =
      sum(
        gene_tpm,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

row_order_df <- bind_rows(
  phylum_summary_order %>%
    mutate(
      order_block = 0
    ),
  genus_order %>%
    mutate(
      order_block = 1
    )
) %>%
  mutate(
    phylum =
      factor(
        phylum,
        levels =
          phylum_order
      )
  ) %>%
  arrange(
    phylum,
    order_block,
    desc(sort_signal),
    desc(total_gene_tpm),
    row_label
  )

row_levels <- row_order_df$row_label

# =========================================================
# SAVE TABLES
# =========================================================

write.table(
  hits_long,
  file.path(
    out_dir,
    "kegg_denitrification_heatmap_orf_hits_with_confidence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  agg,
  file.path(
    out_dir,
    "kegg_denitrification_heatmap_cell_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  rows_kept,
  file.path(
    out_dir,
    "kegg_denitrification_rows_kept.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    phylum = phylum_order
  ),
  file.path(
    out_dir,
    "kegg_denitrification_phylum_order_used.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  row_order_df,
  file.path(
    out_dir,
    "kegg_denitrification_row_order_used.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# MATRIX BUILDERS
# =========================================================

make_abundance_matrix <- function() {

  mat_df <- expand.grid(
    row_label = row_levels,
    gene = gene_order,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      agg %>%
        select(
          row_label,
          gene,
          norm_value
        ),
      by = c(
        "row_label",
        "gene"
      )
    ) %>%
    mutate(
      norm_value =
        replace_na(
          norm_value,
          0
        )
    )

  mat <- mat_df %>%
    select(
      row_label,
      gene,
      norm_value
    ) %>%
    pivot_wider(
      names_from = gene,
      values_from = norm_value
    ) %>%
    as.data.frame()

  rownames(mat) <- mat$row_label
  mat$row_label <- NULL

  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  mat <- mat[
    row_levels,
    gene_order,
    drop = FALSE
  ]

  mat_plot <-
    if (use_log10_for_colour) {
      log10(
        mat +
          plot_pseudocount
      )
    } else {
      mat
    }

  list(
    raw = mat,
    plot = mat_plot
  )
}

make_confidence_label_matrix <- function() {

  conf_df <- expand.grid(
    row_label = row_levels,
    gene = gene_order,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      agg %>%
        mutate(
          confidence_label = ifelse(
            norm_value > 0,
            confidence_label,
            ""
          )
        ) %>%
        select(
          row_label,
          gene,
          confidence_label
        ),
      by = c(
        "row_label",
        "gene"
      )
    ) %>%
    mutate(
      confidence_label =
        replace_na(
          confidence_label,
          ""
        )
    )

  label_mat <- conf_df %>%
    select(
      row_label,
      gene,
      confidence_label
    ) %>%
    pivot_wider(
      names_from = gene,
      values_from =
        confidence_label
    ) %>%
    as.data.frame()

  rownames(label_mat) <-
    label_mat$row_label

  label_mat$row_label <- NULL

  label_mat <- as.matrix(
    label_mat
  )

  label_mat <- label_mat[
    row_levels,
    gene_order,
    drop = FALSE
  ]

  label_mat
}

m_abund <- make_abundance_matrix()
m_label <- make_confidence_label_matrix()

write.table(
  cbind(
    row_label =
      rownames(
        m_abund$raw
      ),
    as.data.frame(
      m_abund$raw
    )
  ),
  file.path(
    out_dir,
    "kegg_denitrification_matrix_abundance_raw.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cbind(
    row_label =
      rownames(
        m_label
      ),
    as.data.frame(
      m_label
    )
  ),
  file.path(
    out_dir,
    "kegg_denitrification_matrix_confidence_letters.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# ANNOTATIONS AND COLOURS
# =========================================================

row_phylum <- row_order_df %>%
  mutate(
    phylum =
      factor(
        phylum,
        levels =
          phylum_order
      )
  ) %>%
  arrange(
    match(
      row_label,
      row_levels
    )
  ) %>%
  pull(phylum)

row_type_vec <- row_order_df %>%
  arrange(
    match(
      row_label,
      row_levels
    )
  ) %>%
  pull(row_type)

row_type_cols <- c(
  "phylum_summary" =
    "#222222",
  "genus" =
    "#4daf4a"
)

left_anno <- rowAnnotation(
  Row_type = row_type_vec,
  col = list(
    Row_type =
      row_type_cols
  ),
  width =
    unit(
      6,
      "mm"
    ),
  show_annotation_name = FALSE,
  annotation_name_gp =
    gpar(
      fontsize = 9,
      fontface = "bold"
    ),
  annotation_legend_param =
    list(
      title = "Row type",
      at =
        names(
          row_type_cols
        ),
      labels =
        c(
          "Phylum",
          "Genus"
        ),
      title_gp =
        gpar(
          fontsize = 9,
          fontface = "bold"
        ),
      labels_gp =
        gpar(
          fontsize = 8
        )
    )
)

finite_vals <- m_abund$plot[
  is.finite(
    m_abund$plot
  )
]

vmin <- min(
  finite_vals,
  na.rm = TRUE
)

vmax <- max(
  finite_vals,
  na.rm = TRUE
)

vmid_candidates <- finite_vals[
  finite_vals > vmin
]

vmid <-
  if (
    length(
      vmid_candidates
    ) > 0
  ) {
    median(
      vmid_candidates,
      na.rm = TRUE
    )
  } else {
    (
      vmin +
        vmax
    ) / 2
  }

abundance_col_fun <- colorRamp2(
  c(
    vmin,
    vmid,
    vmax
  ),
  c(
    "white",
    "gold",
    "firebrick"
  )
)

abundance_legend_title <-
  if (
    use_log10_for_colour
  ) {
    "Relative KEGG gene abundance"
  } else {
    "Normalised TPM"
  }

abundance_legend_breaks <-
  if (
    use_log10_for_colour
  ) {
    log10(
      c(
        1e-6,
        1e-5,
        1e-4,
        1e-3,
        1e-2,
        1e-1
      ) +
        plot_pseudocount
    )
  } else {
    NULL
  }

abundance_legend_labels <-
  if (
    use_log10_for_colour
  ) {
    c(
      "1e-6",
      "1e-5",
      "1e-4",
      "1e-3",
      "1e-2",
      "1e-1"
    )
  } else {
    NULL
  }

col_reactions <- factor(
  gene_to_reaction[
    colnames(
      m_abund$plot
    )
  ],
  levels =
    reaction_levels
)

confidence_legend <- Legend(
  title =
    "KEGG match confidence",
  labels =
    c(
      "H = High",
      "M = Medium",
      "L = Low"
    ),
  type = "points",
  pch =
    c(
      "H",
      "M",
      "L"
    ),
  legend_gp =
    gpar(
      col = "black"
    ),
  labels_gp =
    gpar(
      fontsize = 8
    ),
  title_gp =
    gpar(
      fontsize = 9,
      fontface = "bold"
    )
)

# =========================================================
# HEATMAP
# Colour = taxon-normalised KEGG abundance
# Letter = TPM-weighted average KEGG match confidence
# =========================================================

ht_abundance <- Heatmap(
  m_abund$plot,

  name =
    "KEGG abundance",

  col =
    abundance_col_fun,

  width =
    unit(
      165,
      "mm"
    ),

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_dend = FALSE,
  show_column_dend = FALSE,

  left_annotation =
    left_anno,

  row_split =
    row_phylum,

  row_title = NULL,

  column_split =
    col_reactions,

  row_title_rot = 0,
  row_title_side = "left",

  row_title_gp =
    gpar(
      fontsize = 10,
      fontface = "bold"
    ),

  row_names_side = "left",

  row_names_gp =
    gpar(
      fontsize = 8,
      fontface = ifelse(row_type_vec == "genus", "italic", "plain")
    ),

  row_names_max_width =
    unit(
      65,
      "mm"
    ),

  column_names_gp =
    gpar(
      fontsize = 8,
      fontface = "bold"
    ),

  column_names_rot = 45,
  column_names_side = "bottom",

  column_gap =
    unit(
      4,
      "mm"
    ),

  row_gap =
    unit(
      2,
      "mm"
    ),

  border = TRUE,

  rect_gp =
    gpar(
      col = "grey85",
      lwd = 0.5
    ),

  heatmap_legend_param =
    list(
      title =
        abundance_legend_title,

      at =
        abundance_legend_breaks,

      labels =
        abundance_legend_labels,

      title_gp =
        gpar(
          fontsize = 9,
          fontface = "bold"
        ),

      labels_gp =
        gpar(
          fontsize = 8
        )
    ),

  column_title = expression(
    NO[3]^"-" %->% NO[2]^"-",
    NO[2]^"-" %->% NO,
    NO %->% N[2]*O,
    N[2]*O %->% N[2]
  ),
  column_title_side = "bottom",
  column_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  column_title_rot = 0,

  cell_fun =
    function(
      j,
      i,
      x,
      y,
      width,
      height,
      fill
    ) {

      lab <- m_label[
        i,
        j
      ]

      if (
        !is.na(lab) &&
        lab != ""
      ) {
        grid.text(
          lab,
          x,
          y,
          gp =
            gpar(
              fontsize = 7,
              fontface = "bold",
              col = "black"
            )
        )
      }
    }
)

# =========================================================
# SAVE OUTPUTS
# =========================================================

pdf(
  file.path(
    out_dir,
    "kegg_denitrification_abundance_with_confidence_letters.pdf"
  ),
  width = 19,
  height = 14,
  useDingbats = FALSE
)

draw(
  ht_abundance,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list =
    list(
      confidence_legend
    ),
  padding =
    unit(
      c(
        4,
        4,
        4,
        4
      ),
      "mm"
    )
)

dev.off()

png(
  file.path(
    out_dir,
    "kegg_denitrification_abundance_with_confidence_letters.png"
  ),
  width =
    png_width,
  height =
    png_height,
  res =
    png_res
)

draw(
  ht_abundance,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list =
    list(
      confidence_legend
    ),
  padding =
    unit(
      c(
        4,
        4,
        4,
        4
      ),
      "mm"
    )
)

dev.off()

message("Done.")
message("Output written to: ", out_dir)

message(
  "Main figure: ",
  file.path(
    out_dir,
    "kegg_denitrification_abundance_with_confidence_letters.png"
  )
)

message(
  "ORF validation table: ",
  file.path(
    out_dir,
    "kegg_denitrification_orf_validation_table.tsv"
  )
)

message(
  "Cell confidence table: ",
  file.path(
    out_dir,
    "kegg_denitrification_heatmap_cell_table.tsv"
  )
)
