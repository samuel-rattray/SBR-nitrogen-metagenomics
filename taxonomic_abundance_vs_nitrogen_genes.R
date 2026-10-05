#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SQMtools)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# =========================================================
# USER SETTINGS
# =========================================================

proj <- "/hpcfs/home/w20012148/SqueezeMeta/1st_run"

validation_dir <- file.path(
  proj,
  "results", "paper_results", "fig_3", "now_with_closness", "results"
)

validation_files <- c(
  Nitrification  = "kegg_nitrification_orf_validation_table.tsv",
  Denitrification = "kegg_denitrification_orf_validation_table.tsv",
  DNRA           = "kegg_dnra_orf_validation_table.tsv",
  ANR            = "kegg_assimilatory_nitrate_orf_validation_table.tsv"
)

ranks_to_run <- c("phylum", "class", "order", "family", "genus", "species")

# These are technical DNA-extraction replicates of one sludge sample.
replicate_samples <- c("Sample21", "Sample22", "Sample23")
samples_to_plot <- c(replicate_samples, "mean_all_samples")

# Removed from community tables because they are not taxa.
drop_comm_labels <- c("unmapped", "unassigned", "no hit", "nohit")

# Used only for the optional hidden-unclassified display mode.
hide_taxa_if_contains <- c("unclassified")

out_root <- "scatter_validated_nitrogen_by_rank"

# =========================================================
# COMMAND-LINE INPUT
# Example: Rscript scatter_validated_nitrogen_markers.R 20
# =========================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Provide the number of taxa to plot, e.g. Rscript scatter_validated_nitrogen_markers.R 20")
}

Ntop <- suppressWarnings(as.integer(args[1]))
if (is.na(Ntop) || Ntop < 1) {
  stop("The first argument must be a positive integer.")
}
Nlabel <- Ntop

# =========================================================
# HELPERS
# =========================================================

clean_taxon_names <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unclassified"
  trimws(x)
}

norm_taxon <- function(x) tolower(trimws(as.character(x)))

is_hidden_taxon <- function(x, patterns = character(0)) {
  if (length(patterns) == 0) return(rep(FALSE, length(x)))
  grepl(paste(patterns, collapse = "|"), norm_taxon(x), ignore.case = TRUE)
}

resolve_abund_col <- function(sample_name, cols) {
  if (sample_name %in% cols) return(sample_name)

  hits <- grep(sample_name, cols, value = TRUE, fixed = TRUE)
  if (length(hits) == 0) {
    stop(
      "Could not find an ORF abundance column for '", sample_name,
      "'. Available columns: ", paste(cols, collapse = ", ")
    )
  }

  tpm <- hits[grepl("TPM", hits, ignore.case = TRUE)]
  if (length(tpm) > 0) return(tpm[1])

  cov <- hits[grepl("Coverage", hits, ignore.case = TRUE)]
  if (length(cov) > 0) return(cov[1])

  rpkm <- hits[grepl("RPKM", hits, ignore.case = TRUE)]
  if (length(rpkm) > 0) return(rpkm[1])

  hits[1]
}

# =========================================================
# LOAD SQM PROJECT
# =========================================================

SQM <- loadSQM(proj, tax_mode = "prokfilter")

orf_tax   <- SQM$orfs$tax
orf_annot <- SQM$orfs$table
orf_abund <- SQM$orfs$abund

if (is.null(orf_tax)) stop("SQM$orfs$tax is NULL.")
if (is.null(orf_annot)) stop("SQM$orfs$table is NULL.")
if (is.null(orf_abund)) stop("SQM$orfs$abund is NULL.")

all_orf_ids_global <- Reduce(
  intersect,
  list(rownames(orf_tax), rownames(orf_annot), rownames(orf_abund))
)
if (length(all_orf_ids_global) == 0) {
  stop("No overlapping ORF IDs between SQM taxonomy, annotation and abundance tables.")
}

resolved_sample_cols <- setNames(
  vapply(
    replicate_samples,
    resolve_abund_col,
    FUN.VALUE = character(1),
    cols = colnames(orf_abund)
  ),
  replicate_samples
)

message(
  "ORF abundance columns used: ",
  paste(paste(names(resolved_sample_cols), resolved_sample_cols, sep = "="), collapse = "; ")
)

# =========================================================
# LOAD VALIDATED NITROGEN-MARKER ORFS
#
# This replaces the old broad KEGGPATH == 'Nitrogen metabolism' selection.
# Only H/M/L calls from the four pathway validation scripts are retained.
# U calls are excluded. If an ORF occurs in multiple pathway tables, it is
# counted once in this broad scatter-plot numerator.
# =========================================================

validation_list <- list()

for (pathway_name in names(validation_files)) {
  f <- file.path(validation_dir, validation_files[[pathway_name]])

  if (!file.exists(f)) {
    stop(
      "Missing validation table: ", f, "\n",
      "Rerun the corresponding final pathway heatmap script first."
    )
  }

  x <- fread(f, data.table = FALSE)
  needed <- c("orf_id", "gene", "confidence")
  missing <- setdiff(needed, colnames(x))

  if (length(missing) > 0) {
    stop(
      basename(f), " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }

  x$pathway <- pathway_name
  validation_list[[pathway_name]] <- x
}

validation_all <- rbindlist(validation_list, fill = TRUE, use.names = TRUE)
validation_all$orf_id <- as.character(validation_all$orf_id)
validation_all$gene <- as.character(validation_all$gene)
validation_all$confidence <- as.character(validation_all$confidence)
validation_all$pathway <- as.character(validation_all$pathway)

validated_calls <- validation_all[
  validation_all$confidence %in% c("High", "Medium", "Low"),
  , drop = FALSE
]

if (nrow(validated_calls) == 0) {
  stop("No H/M/L validated nitrogen-marker ORFs were found.")
}

# Keep one row per ORF for counting, but retain provenance in an output table.
conf_rank <- c("Low" = 1, "Medium" = 2, "High" = 3)

validated_catalog <- validated_calls %>%
  group_by(orf_id) %>%
  summarise(
    pathways = paste(sort(unique(pathway)), collapse = ";"),
    genes = paste(sort(unique(gene)), collapse = ";"),
    best_confidence = names(conf_rank)[
      which.max(vapply(names(conf_rank), function(z) any(confidence == z), logical(1)) * conf_rank)
    ],
    n_validated_calls = n(),
    .groups = "drop"
  )

validated_orf_ids_all <- unique(validated_catalog$orf_id)
validated_orf_ids_global <- intersect(validated_orf_ids_all, all_orf_ids_global)
missing_validated_ids <- setdiff(validated_orf_ids_all, all_orf_ids_global)

message("Unique H/M/L validated marker ORFs: ", length(validated_orf_ids_all))
message("Validated marker ORFs found in SQM project: ", length(validated_orf_ids_global))
message("Validated marker ORFs absent from SQM project: ", length(missing_validated_ids))

if (length(validated_orf_ids_global) == 0) {
  stop("None of the validated marker ORF IDs matched the SQM project ORF IDs.")
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write.table(
  validated_catalog,
  file.path(out_root, "validated_nitrogen_marker_catalog.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (length(missing_validated_ids) > 0) {
  write.table(
    data.frame(orf_id = missing_validated_ids),
    file.path(out_root, "validated_ORFs_not_found_in_SQM.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# =========================================================
# TAXONOMY SETUP
# =========================================================

phylum_candidates <- c("phylum", "Phylum", "p", "tax_phylum")
phylum_col <- intersect(phylum_candidates, colnames(orf_tax))
if (length(phylum_col) == 0) {
  stop("Could not find a phylum column in SQM$orfs$tax.")
}
phylum_col <- phylum_col[1]

all_phyla_global <- sort(unique(clean_taxon_names(
  orf_tax[all_orf_ids_global, phylum_col]
)))

phylum_palette <- grDevices::hcl.colors(length(all_phyla_global), palette = "Dark 3")
phylum_colour_map <- stats::setNames(phylum_palette, all_phyla_global)
unknown_phylum_colour <- "grey70"

# =========================================================
# RUN ONE TAXONOMIC RANK
# =========================================================

run_one_rank <- function(rank_name) {

  message("========================================")
  message("Running rank: ", rank_name)

  rank_candidates <- c(
    rank_name,
    tools::toTitleCase(rank_name),
    substr(rank_name, 1, 1),
    paste0("tax_", rank_name)
  )

  rank_col <- intersect(rank_candidates, colnames(orf_tax))
  if (length(rank_col) == 0) {
    warning("Skipping rank '", rank_name, "': taxonomy column not found.")
    return(NULL)
  }
  rank_col <- rank_col[1]

  if (
    is.null(SQM$taxa) ||
    is.null(SQM$taxa[[rank_name]]) ||
    is.null(SQM$taxa[[rank_name]]$percent)
  ) {
    warning("Skipping rank '", rank_name, "': community percent table not found.")
    return(NULL)
  }

  comm_mat <- SQM$taxa[[rank_name]]$percent
  available_comm_samples <- colnames(comm_mat)

  missing_comm <- setdiff(replicate_samples, available_comm_samples)
  if (length(missing_comm) > 0) {
    stop("Missing community sample columns: ", paste(missing_comm, collapse = ", "))
  }

  build_one_sample_table <- function(sample_name) {

    # -------------------------
    # Community abundance
    # -------------------------
    if (sample_name == "mean_all_samples") {
      community_vals <- rowMeans(
        comm_mat[, replicate_samples, drop = FALSE],
        na.rm = TRUE
      )
    } else {
      community_vals <- comm_mat[, sample_name]
    }

    community <- data.frame(
      taxon = clean_taxon_names(rownames(comm_mat)),
      community_percent = as.numeric(community_vals),
      stringsAsFactors = FALSE
    )
    community$community_percent[!is.finite(community$community_percent)] <- 0
    community <- community[
      !(norm_taxon(community$taxon) %in% norm_taxon(drop_comm_labels)),
      , drop = FALSE
    ]
    community <- aggregate(community_percent ~ taxon, community, sum)

    # -------------------------
    # ORF abundance
    # mean_all_samples now uses ONLY Sample21-23, not every column in SQM$orfs$abund.
    # -------------------------
    if (sample_name == "mean_all_samples") {
      all_abund <- rowMeans(
        orf_abund[all_orf_ids_global, unname(resolved_sample_cols), drop = FALSE],
        na.rm = TRUE
      )
      nit_abund <- rowMeans(
        orf_abund[validated_orf_ids_global, unname(resolved_sample_cols), drop = FALSE],
        na.rm = TRUE
      )
      abundance_column_used <- paste(unname(resolved_sample_cols), collapse = ";")
    } else {
      ac <- resolved_sample_cols[[sample_name]]
      all_abund <- as.numeric(orf_abund[all_orf_ids_global, ac])
      nit_abund <- as.numeric(orf_abund[validated_orf_ids_global, ac])
      abundance_column_used <- ac
    }

    # -------------------------
    # All ORFs by taxon
    # -------------------------
    df_all <- data.frame(
      orf = all_orf_ids_global,
      abund = as.numeric(all_abund),
      taxon = clean_taxon_names(orf_tax[all_orf_ids_global, rank_col]),
      phylum = clean_taxon_names(orf_tax[all_orf_ids_global, phylum_col]),
      stringsAsFactors = FALSE
    )
    df_all <- df_all[is.finite(df_all$abund) & df_all$abund > 0, , drop = FALSE]

    if (nrow(df_all) == 0) {
      stop("No positive-abundance ORFs for ", sample_name, " at rank ", rank_name)
    }

    all_taxon_abs <- aggregate(abund ~ taxon, df_all, sum)
    colnames(all_taxon_abs)[2] <- "all_orf_abund"

    # Use the abundance-dominant phylum when the same taxon label maps to >1 phylum.
    taxon_phylum_abs <- aggregate(abund ~ taxon + phylum, df_all, sum)
    taxon_phylum_abs <- taxon_phylum_abs[
      order(taxon_phylum_abs$taxon, -taxon_phylum_abs$abund, taxon_phylum_abs$phylum),
    ]
    taxon_to_phylum <- taxon_phylum_abs[
      !duplicated(taxon_phylum_abs$taxon), c("taxon", "phylum")
    ]
    if (rank_name == "phylum") taxon_to_phylum$phylum <- taxon_to_phylum$taxon

    # -------------------------
    # Validated H/M/L nitrogen-marker ORFs by taxon
    # -------------------------
    df_nit <- data.frame(
      orf = validated_orf_ids_global,
      abund = as.numeric(nit_abund),
      taxon = clean_taxon_names(orf_tax[validated_orf_ids_global, rank_col]),
      phylum = clean_taxon_names(orf_tax[validated_orf_ids_global, phylum_col]),
      stringsAsFactors = FALSE
    )
    df_nit <- df_nit[is.finite(df_nit$abund) & df_nit$abund > 0, , drop = FALSE]

    if (nrow(df_nit) > 0) {
      nit_taxon_abs <- aggregate(abund ~ taxon, df_nit, sum)
      colnames(nit_taxon_abs)[2] <- "validated_nitrogen_orf_abund"
    } else {
      nit_taxon_abs <- data.frame(
        taxon = character(0),
        validated_nitrogen_orf_abund = numeric(0)
      )
    }

    # -------------------------
    # Merge and normalise within taxon
    # -------------------------
    merged <- merge(community, all_taxon_abs, by = "taxon", all = TRUE)
    merged <- merge(merged, nit_taxon_abs, by = "taxon", all = TRUE)
    merged <- merge(merged, taxon_to_phylum, by = "taxon", all.x = TRUE)

    merged$community_percent[is.na(merged$community_percent)] <- 0
    merged$all_orf_abund[is.na(merged$all_orf_abund)] <- 0
    merged$validated_nitrogen_orf_abund[is.na(merged$validated_nitrogen_orf_abund)] <- 0

    if (rank_name == "phylum") merged$phylum <- merged$taxon

    merged$validated_nitrogen_within_taxon_percent <- ifelse(
      merged$all_orf_abund > 0,
      100 * merged$validated_nitrogen_orf_abund / merged$all_orf_abund,
      0
    )

    total_all <- sum(merged$all_orf_abund, na.rm = TRUE)
    total_nit <- sum(merged$validated_nitrogen_orf_abund, na.rm = TRUE)
    overall_nit_percent <- ifelse(total_all > 0, 100 * total_nit / total_all, NA_real_)

    merged$validated_nitrogen_enrichment_ratio <- ifelse(
      is.finite(overall_nit_percent) & overall_nit_percent > 0,
      merged$validated_nitrogen_within_taxon_percent / overall_nit_percent,
      NA_real_
    )

    merged$sample <- sample_name
    merged$abund_column_used <- abundance_column_used
    merged$is_hidden_unclassified <- is_hidden_taxon(
      merged$taxon, hide_taxa_if_contains
    )

    merged
  }

  # -------------------------
  # Build the three technical-replicate tables.
  # No inferential tests are performed because they are technical replicates.
  # -------------------------
  replicate_tables <- lapply(replicate_samples, build_one_sample_table)
  names(replicate_tables) <- replicate_samples
  replicate_long_all <- do.call(rbind, replicate_tables)
  rownames(replicate_long_all) <- NULL

  # Top taxa are selected from mean COMMUNITY abundance only.
  mean_comm_all <- aggregate(
    community_percent ~ taxon,
    replicate_long_all,
    mean,
    na.rm = TRUE
  )

  phylum_lookup <- unique(
    replicate_long_all[, c("taxon", "phylum", "is_hidden_unclassified")]
  )
  phylum_lookup <- phylum_lookup[!duplicated(phylum_lookup$taxon), , drop = FALSE]
  mean_comm_all <- merge(mean_comm_all, phylum_lookup, by = "taxon", all.x = TRUE)
  mean_comm_all <- mean_comm_all[
    order(mean_comm_all$community_percent, decreasing = TRUE), , drop = FALSE
  ]

  # Descriptive mean +/- SD across technical replicates.
  make_desc <- function(taxon_name) {
    z <- replicate_long_all[replicate_long_all$taxon == taxon_name, , drop = FALSE]
    data.frame(
      taxon = taxon_name,
      community_percent_mean = mean(z$community_percent, na.rm = TRUE),
      community_percent_sd = sd(z$community_percent, na.rm = TRUE),
      validated_nitrogen_percent_mean = mean(z$validated_nitrogen_within_taxon_percent, na.rm = TRUE),
      validated_nitrogen_percent_sd = sd(z$validated_nitrogen_within_taxon_percent, na.rm = TRUE),
      validated_nitrogen_orf_abund_mean = mean(z$validated_nitrogen_orf_abund, na.rm = TRUE),
      validated_nitrogen_orf_abund_sd = sd(z$validated_nitrogen_orf_abund, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  descriptive_summary <- do.call(
    rbind,
    lapply(unique(replicate_long_all$taxon), make_desc)
  )
  descriptive_summary <- merge(
    descriptive_summary,
    phylum_lookup[, c("taxon", "phylum", "is_hidden_unclassified")],
    by = "taxon", all.x = TRUE
  )

  modes <- list(
    show_unclassified = FALSE,
    hide_unclassified = TRUE
  )

  for (mode_name in names(modes)) {

    hide_mode <- modes[[mode_name]]
    rank_out_dir <- file.path(out_root, mode_name, rank_name)
    dir.create(rank_out_dir, recursive = TRUE, showWarnings = FALSE)

    mean_comm_mode <- mean_comm_all
    if (hide_mode) {
      mean_comm_mode <- mean_comm_mode[
        !mean_comm_mode$is_hidden_unclassified, , drop = FALSE
      ]
    }

    mean_comm_mode <- mean_comm_mode[
      order(mean_comm_mode$community_percent, decreasing = TRUE), , drop = FALSE
    ]
    top_taxa <- head(mean_comm_mode$taxon, Ntop)
    top_taxa_ordered <- mean_comm_mode$taxon[mean_comm_mode$taxon %in% top_taxa]

    make_one_plot <- function(sample_name) {

      plot_df <- build_one_sample_table(sample_name)
      plot_df <- plot_df[plot_df$taxon %in% top_taxa, , drop = FALSE]

      if (hide_mode) {
        plot_df <- plot_df[!plot_df$is_hidden_unclassified, , drop = FALSE]
      }
      if (nrow(plot_df) == 0) return(NULL)

      plot_df$taxon <- factor(as.character(plot_df$taxon), levels = top_taxa_ordered)
      plot_df <- plot_df[order(plot_df$taxon), , drop = FALSE]
      plot_df$label <- as.character(plot_df$taxon)

      phyla_in_plot <- sort(unique(na.omit(as.character(plot_df$phylum))))
      x_max <- max(plot_df$community_percent, na.rm = TRUE)
      y_max <- max(plot_df$validated_nitrogen_within_taxon_percent, na.rm = TRUE)
      if (!is.finite(x_max) || x_max <= 0) x_max <- 1
      if (!is.finite(y_max) || y_max <= 0) y_max <- 1

      subtitle_text <- paste0(
        "Top ", nrow(plot_df), " ", rank_name,
        " taxa by mean community abundance across technical extraction replicates; ",
        "point size = nitrogen gene abundance"
      )
      if (hide_mode) {
        subtitle_text <- paste0(
          subtitle_text,
          "; unclassified taxa retained in calculations but hidden from display"
        )
      }

      p <- ggplot(
        plot_df,
        aes(
          x = community_percent,
          y = validated_nitrogen_within_taxon_percent,
          size = validated_nitrogen_orf_abund
        )
      ) +
        geom_point(
          aes(fill = phylum),
          shape = 21,
          colour = "black",
          alpha = 0.9,
          stroke = 0.5
        ) +
        geom_text_repel(
          aes(label = label),
          size = 3,
          box.padding = 0.35,
          point.padding = 0.25,
          force = 1,
          force_pull = 0.5,
          max.iter = 10000,
          max.overlaps = Inf,
          min.segment.length = 0,
          segment.alpha = 0.6,
          seed = 123
        ) +
        scale_fill_manual(
          values = phylum_colour_map,
          breaks = phyla_in_plot,
          drop = FALSE,
          na.value = unknown_phylum_colour
        ) +
        theme_bw(base_size = 11) +
        labs(
          title = paste(
            "Community abundance vs validated nitrogen-transformation markers:",
            sample_name
          ),
          subtitle = subtitle_text,
          x = paste0(rank_name, " abundance in community (%)"),
          y = paste0("Nitrogen metabolism gene abundance within ", rank_name, " (%)"),
          size = "Nitrogen gene\nabundance",
          fill = "Phylum"
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold"),
          plot.margin = margin(t = 5.5, r = 130, b = 5.5, l = 5.5)
        ) +
        coord_cartesian(
          xlim = c(0, x_max * 1.55),
          ylim = c(0, y_max * 1.12),
          clip = "off"
        )

      suffix_tag <- if (hide_mode) "hidden_unclassified" else "shown_unclassified"
      out_prefix <- file.path(
        rank_out_dir,
        paste0(
          "FigureX_Scatter_Community_vs_ValidatedNitrogen_",
          rank_name, "_top", Ntop, "_", suffix_tag, "_", sample_name
        )
      )

      ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 6.5, dpi = 300)
      ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = 6.5)

      write.table(
        plot_df,
        paste0(out_prefix, "_data.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }

    for (s in samples_to_plot) make_one_plot(s)

    # Descriptive outputs only; no Friedman/Wilcoxon testing.
    write.table(
      replicate_long_all,
      file.path(rank_out_dir, paste0("ValidatedNitrogen_", rank_name, "_all_taxa_technical_replicates.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
      descriptive_summary[descriptive_summary$taxon %in% top_taxa, , drop = FALSE],
      file.path(rank_out_dir, paste0("ValidatedNitrogen_", rank_name, "_top", Ntop, "_mean_SD.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
      data.frame(taxon = top_taxa_ordered),
      file.path(rank_out_dir, paste0("ValidatedNitrogen_", rank_name, "_top", Ntop, "_selected_taxa.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    message("ALL DONE FOR RANK: ", rank_name, " | MODE: ", mode_name)
  }

  invisible(TRUE)
}

# =========================================================
# RUN
# =========================================================

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (this_rank in ranks_to_run) {
  tryCatch(
    run_one_rank(this_rank),
    error = function(e) {
      message("ERROR while running rank ", this_rank, ": ", conditionMessage(e))
    }
  )
}

message("========================================")
message("FINISHED ALL REQUESTED RANKS")
message("Output root: ", out_root)

