#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sgtr)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
  library(scales)
  library(png)
  library(grid)
})

# -----------------------------
# Helpers
# -----------------------------
die <- function(msg) {
  cat("\nERROR:", msg, "\n", file = stderr())
  quit(status = 1)
}

info <- function(msg) cat(msg, "\n")

require_file <- function(fname) {
  if (!file.exists(fname)) {
    die(paste0(
      "This folder must contain the file '", fname, "'.\n",
      "Please run this script inside the PSASS analyze folder (e.g., .../03_analyze/50k/)."
    ))
  }
}

# Parse CLI args: --n 25
get_arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(NA_character_)
  idx <- which(args == flag)
  if (length(idx) == 0) return(NA_character_)
  if (idx[1] == length(args)) return(NA_character_)
  args[idx[1] + 1]
}

parse_n_from_cli <- function() {
  v <- get_arg_value("--n")
  if (is.na(v)) return(NA_integer_)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 1) {
    die("Invalid value for --n. Use a positive integer, e.g. --n 25")
  }
  n
}

# Read a table that may be tab-separated OR whitespace-separated
read_flexible <- function(path) {
  df <- tryCatch(
    read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df2 <- tryCatch(
    read.table(path, sep = "", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df2) && ncol(df2) >= 2) return(df2)

  die(paste0("Could not parse input file: ", path, " (neither tab nor whitespace separated)."))
}

# Remove duplicated header line if present
fix_psass_window <- function(infile, outfile) {
  lines <- readLines(infile, warn = FALSE)
  if (length(lines) < 2) die(paste0("Input looks empty: ", infile))

  if (lines[1] == lines[2]) {
    info("Detected duplicated header line in psass_window.tsv -> removing it.")
    lines <- lines[c(1, 3:length(lines))]
  } else {
    info("No duplicated header detected in psass_window.tsv.")
  }

  writeLines(lines, outfile)
  info(paste0("Wrote cleaned window file: ", outfile))
}

# Detect prefix for outputs based on folder path
detect_prefix <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (grepl("psass_hap1", wd, ignore.case = TRUE)) return("HAP1_PSASS")
  if (grepl("psass_hap2", wd, ignore.case = TRUE)) return("HAP2_PSASS")
  return("PSASS")
}

# Order contigs in an "obvious" way
order_contigs_obvious <- function(contig_sizes_df) {
  contigs <- contig_sizes_df$Contig

  chr_num <- suppressWarnings(as.integer(sub("^Chromosome", "", contigs)))
  if (length(contigs) > 0 && all(!is.na(chr_num))) {
    return(contigs[order(chr_num)])
  }

  sca_num <- suppressWarnings(as.integer(sub("^(scaffold|Scaffold)[_\\-]?", "", contigs)))
  if (length(contigs) > 0 && all(!is.na(sca_num))) {
    return(contigs[order(sca_num)])
  }

  contig_sizes_df <- contig_sizes_df[order(contig_sizes_df$Length, decreasing = TRUE), ]
  contig_sizes_df$Contig
}

# Write chromosomes.tsv mapping Contig -> Label
write_chromosomes_tsv <- function(contigs_ordered, outfile) {
  labels <- contigs_ordered
  labels <- sub("^Chromosome", "Chr", labels)
  labels <- sub("^(scaffold|Scaffold)[_\\-]?", "Scf", labels)

  out <- data.frame(Contig = contigs_ordered, Label = labels, stringsAsFactors = FALSE)
  write.table(out, outfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  info(paste0("Wrote chromosomes mapping file: ", outfile))
}

# -----------------------------
# Add center legend to circos PNG
# -----------------------------
add_circos_center_legend_png <- function(input_png, output_png) {
  img <- png::readPNG(input_png)
  h <- dim(img)[1]
  w <- dim(img)[2]

  grDevices::png(output_png, width = w, height = h, res = 150)
  grid::grid.newpage()
  grid::grid.raster(img, x = 0.5, y = 0.5, width = 1, height = 1)

  # caixa central
  box_w <- 0.22
  box_h <- 0.16
  box_x <- 0.50
  box_y <- 0.50

  grid::grid.roundrect(
    x = box_x, y = box_y,
    width = box_w, height = box_h,
    r = unit(0.02, "snpc"),
    gp = grid::gpar(
      fill = rgb(1, 1, 1, 0.82),
      col = "grey65",
      lwd = 1.2
    )
  )

  grid::grid.text(
    "Legend",
    x = box_x, y = box_y + 0.055,
    gp = grid::gpar(fontsize = 15, fontface = "bold", col = "black")
  )

  legend_y <- c(0.52, 0.495, 0.47, 0.445)
  legend_labels <- c("Fst", "Female SNPs", "Male SNPs", "Depth ratio")
  legend_cols <- c(COL_GREEN_DARK, COL_PINK_DARK, COL_BLUE_DARK, COL_PURPLE_DARK)

  for (i in seq_along(legend_y)) {
    grid::grid.rect(
      x = 0.445, y = legend_y[i],
      width = 0.012, height = 0.012,
      gp = grid::gpar(fill = legend_cols[i], col = legend_cols[i])
    )
    grid::grid.text(
      legend_labels[i],
      x = 0.458, y = legend_y[i],
      just = "left",
      gp = grid::gpar(fontsize = 11, col = "black")
    )
  }

  grDevices::dev.off()
  info(paste0("Wrote circos with center legend: ", output_png))
}

# -----------------------------
# Manhattan helpers
# -----------------------------
build_manhattan_plot <- function(df_filtered, keep_contigs, shared_snps_max, fst_raw_max, fst_gt1_n, prefix) {
  df_plot <- df_filtered %>%
    mutate(
      Position_Mbp = Position / 1e6,
      Depth_ratio_plot = pmin(Depth_ratio, 5)
    )

  chrom_lengths <- df_plot %>%
    group_by(Chromosome = Contig) %>%
    summarise(
      Length_bp = max(Position, na.rm = TRUE),
      Length_Mbp = max(Position, na.rm = TRUE) / 1e6,
      .groups = "drop"
    ) %>%
    arrange(factor(Chromosome, levels = keep_contigs)) %>%
    mutate(
      Cumulative_start_bp = cumsum(dplyr::lag(Length_bp, default = 0)),
      Cumulative_start_Mbp = Cumulative_start_bp / 1e6,
      Cumulative_end_Mbp = Cumulative_start_Mbp + Length_Mbp,
      Midpoint_Mbp = (Cumulative_start_Mbp + Cumulative_end_Mbp) / 2,
      Color_index = ifelse(row_number() %% 2 == 0, "even", "odd")
    )

  df_plot <- df_plot %>%
    left_join(
      chrom_lengths %>% select(Chromosome, Cumulative_start_Mbp, Color_index),
      by = c("Contig" = "Chromosome")
    ) %>%
    mutate(
      Cumulative_position_Mbp = Cumulative_start_Mbp + Position_Mbp
    )

  y_max_snps <- shared_snps_max * 1.05
  y_max_fst  <- min(fst_raw_max * 1.05, 1)
  y_max_depth_ratio <- 5

  chr_labels <- gsub("^Chromosome", "Chr", chrom_lengths$Chromosome)
  chr_labels <- gsub("^(scaffold|Scaffold)[_\\-]?", "Scf", chr_labels)

  info("Generating Manhattan plot with ggplot2 ...")
  info(paste0("  Manhattan SNP Y max (shared): ", round(y_max_snps, 2)))
  info(paste0("  Manhattan Fst display max: ", round(y_max_fst, 4)))
  info("  Manhattan Depth_ratio display max: 5")

  create_male_plot <- function(data, y_max) {
    ggplot(data) +
      geom_point(
        aes(x = Cumulative_position_Mbp, y = Snps_males, color = Color_index),
        size = 0.3,
        alpha = 0.8,
        shape = 16
      ) +
      geom_smooth(
        aes(x = Cumulative_position_Mbp, y = Snps_males, group = Contig),
        method = "loess",
        se = FALSE,
        color = COL_BLUE_DARK,
        linewidth = 1.0,
        span = 0.05
      ) +
      scale_color_manual(
        values = c("odd" = COL_BLUE_LIGHT, "even" = COL_BLUE_DARK),
        guide = "none"
      ) +
      scale_x_continuous(
        name = NULL,
        breaks = chrom_lengths$Midpoint_Mbp,
        labels = chr_labels,
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(
        name = "Male SNPs (per window)",
        limits = c(0, y_max),
        breaks = pretty_breaks(n = 5),
        labels = comma,
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      labs(
        title = paste("Male-specific SNPs (Max:",
                      round(max(data$Snps_males, na.rm = TRUE), 0), ")")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(5, 10, 0, 10)
      )
  }

  create_female_plot <- function(data, y_max) {
    ggplot(data) +
      geom_point(
        aes(x = Cumulative_position_Mbp, y = Snps_females, color = Color_index),
        size = 0.3,
        alpha = 0.8,
        shape = 16
      ) +
      geom_smooth(
        aes(x = Cumulative_position_Mbp, y = Snps_females, group = Contig),
        method = "loess",
        se = FALSE,
        color = COL_PINK_DARK,
        linewidth = 1.0,
        span = 0.05
      ) +
      scale_color_manual(
        values = c("odd" = COL_PINK_LIGHT, "even" = COL_PINK_DARK),
        guide = "none"
      ) +
      scale_x_continuous(
        name = NULL,
        breaks = chrom_lengths$Midpoint_Mbp,
        labels = chr_labels,
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(
        name = "Female SNPs (per window)",
        limits = c(0, y_max),
        breaks = pretty_breaks(n = 5),
        labels = comma,
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      labs(
        title = paste("Female-specific SNPs (Max:",
                      round(max(data$Snps_females, na.rm = TRUE), 0), ")")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(0, 10, 0, 10)
      )
  }

  create_fst_plot <- function(data, y_max) {
    ggplot(data) +
      geom_point(
        aes(x = Cumulative_position_Mbp, y = Fst, color = Color_index),
        size = 0.3,
        alpha = 0.8,
        shape = 16
      ) +
      geom_smooth(
        aes(x = Cumulative_position_Mbp, y = Fst, group = Contig),
        method = "loess",
        se = FALSE,
        color = COL_GREEN_DARK,
        linewidth = 1.0,
        span = 0.05
      ) +
      geom_hline(
        yintercept = 0.25,
        linetype = "dashed",
        color = "gray30",
        alpha = 0.8
      ) +
      scale_color_manual(
        values = c("odd" = COL_GREEN_LIGHT, "even" = COL_GREEN_DARK),
        guide = "none"
      ) +
      scale_x_continuous(
        name = NULL,
        breaks = chrom_lengths$Midpoint_Mbp,
        labels = chr_labels,
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(
        name = "Fst",
        limits = c(0, y_max),
        breaks = pretty_breaks(n = 5),
        labels = number_format(accuracy = 0.01),
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      labs(
        title = paste("Fst (Raw max:", round(fst_raw_max, 2),
                      "| >1:", fst_gt1_n, "windows)")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(0, 10, 0, 10)
      )
  }

  create_depth_ratio_plot <- function(data, y_max) {
    ggplot(data) +
      geom_point(
        aes(x = Cumulative_position_Mbp, y = Depth_ratio_plot, color = Color_index),
        size = 0.3,
        alpha = 0.8,
        shape = 16
      ) +
      geom_smooth(
        aes(x = Cumulative_position_Mbp, y = Depth_ratio_plot, group = Contig),
        method = "loess",
        se = FALSE,
        color = COL_PURPLE_DARK,
        linewidth = 1.0,
        span = 0.05
      ) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "gray30",
        alpha = 0.8
      ) +
      scale_color_manual(
        values = c("odd" = COL_PURPLE_LIGHT, "even" = COL_PURPLE_DARK),
        guide = "none"
      ) +
      scale_x_continuous(
        name = "Chromosome",
        breaks = chrom_lengths$Midpoint_Mbp,
        labels = chr_labels,
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(
        name = "Depth ratio",
        limits = c(0, y_max),
        breaks = pretty_breaks(n = 5),
        labels = number_format(accuracy = 0.01),
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      labs(
        title = paste("Depth ratio (display max:", y_max,
                      "| raw max:", round(max(data$Depth_ratio, na.rm = TRUE), 2), ")")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(0, 10, 10, 10)
      )
  }

  p_male   <- create_male_plot(df_plot, y_max_snps)
  p_female <- create_female_plot(df_plot, y_max_snps)
  p_fst    <- create_fst_plot(df_plot, y_max_fst)
  p_depth  <- create_depth_ratio_plot(df_plot, y_max_depth_ratio)

  combined_plot <- p_male / p_female / p_fst / p_depth +
    plot_layout(heights = c(1, 1, 1.2, 1.1)) +
    plot_annotation(
      title = "PSASS Analysis: SNP counts, Fst distribution and depth ratio",
      subtitle = paste(
        "Scaffolds:", length(keep_contigs),
        "| Same Y scale for SNPs (max:", round(shared_snps_max, 0), ")"
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )

  out_png <- paste0(prefix, "_manhattan_FST_SNPf_SNPm_pastelAlt.png")

  ggsave(
    out_png,
    combined_plot,
    width = 16,
    height = 17,
    device = "png",
    dpi = 300,
    bg = "white"
  )

  info(paste0("Wrote Manhattan plot: ", out_png))
}

# -----------------------------
# Main
# -----------------------------
REQ_WINDOW   <- "psass_window.tsv"
CLEAN_WINDOW <- "psass_window.clean.tsv"
FILTERED_WIN <- "psass_window.filtered.tsv"
CHR_TSV      <- "chromosomes.tsv"
SELECTED_TXT <- "selected_contigs.txt"

REQUESTED_N <- parse_n_from_cli()

info("------------------------------------------------------------")
info("PSASS plotting helper")
info("Input required in this folder: psass_window.tsv")
info("Usage examples:")
info("  Rscript run_psass_plot_auto.R --n 25")
info("  Rscript run_psass_plot_auto.R            (uses ALL contigs)")
info("------------------------------------------------------------")

require_file(REQ_WINDOW)

info("Reading psass_window.tsv ...")
fix_psass_window(REQ_WINDOW, CLEAN_WINDOW)

df <- read_flexible(CLEAN_WINDOW)

needed_cols <- c("Contig", "Position", "Fst", "Snps_females", "Snps_males", "Length", "Depth_ratio")
missing <- setdiff(needed_cols, colnames(df))
if (length(missing) > 0) {
  die(paste0(
    "psass_window.clean.tsv is missing required columns: ",
    paste(missing, collapse = ", "),
    "\nExpected at least: ", paste(needed_cols, collapse = ", ")
  ))
}

# numeric coercion
df$Position      <- suppressWarnings(as.numeric(df$Position))
df$Fst           <- suppressWarnings(as.numeric(df$Fst))
df$Snps_females  <- suppressWarnings(as.numeric(df$Snps_females))
df$Snps_males    <- suppressWarnings(as.numeric(df$Snps_males))
df$Length        <- suppressWarnings(as.numeric(df$Length))
df$Depth_ratio   <- suppressWarnings(as.numeric(df$Depth_ratio))

contig_sizes <- aggregate(Length ~ Contig, data = df, FUN = max)
ordered_contigs <- order_contigs_obvious(contig_sizes)
max_n <- length(ordered_contigs)

if (is.na(REQUESTED_N)) {
  n_keep <- max_n
  info(paste0("No --n provided. Keeping ALL contigs for analysis (", max_n, ")."))
  info("Tip: if your assembly is fragmented (many scaffolds), use e.g. --n 25 or --n 50 to avoid Circos errors.")
} else {
  n_keep <- min(REQUESTED_N, max_n)
  if (REQUESTED_N > max_n) {
    info(paste0("You requested ", REQUESTED_N, " contigs, but only ", max_n, " are available. Using ALL (", max_n, ")."))
  } else {
    info(paste0("Keeping the first ", n_keep, " contigs for analysis."))
  }
}

keep_contigs <- head(ordered_contigs, n_keep)

writeLines(keep_contigs, SELECTED_TXT)
info(paste0("Wrote selected contigs list: ", SELECTED_TXT))

df_filt <- df[df$Contig %in% keep_contigs, ]

if (nrow(df_filt) == 0) {
  die("After filtering contigs, no rows remained in the dataset.")
}

# Fst do circos truncado entre 0 e 1 para evitar picos estranhos e saídas do cromossomo
df_filt$Fst_circos <- pmax(pmin(df_filt$Fst, 1), 0)

# Depth ratio truncado para visualização
df_filt$Depth_ratio_circos <- pmin(df_filt$Depth_ratio, 5)

write.table(df_filt, FILTERED_WIN, sep = "\t", row.names = FALSE, quote = FALSE)
info(paste0("Wrote filtered window file: ", FILTERED_WIN))

write_chromosomes_tsv(keep_contigs, CHR_TSV)

# -----------------------------
# Scales
# -----------------------------
shared_snps_max <- max(c(df_filt$Snps_females, df_filt$Snps_males), na.rm = TRUE)

if (!is.finite(shared_snps_max) || is.na(shared_snps_max) || shared_snps_max <= 0) {
  die("Could not determine a valid shared SNP Y limit from Snps_females / Snps_males.")
}

fst_raw_max <- max(df_filt$Fst, na.rm = TRUE)
fst_gt1_n   <- sum(df_filt$Fst > 1, na.rm = TRUE)

info(paste0("Shared SNP y-limit: 0 to ", shared_snps_max))
info(paste0("Raw Fst max in file: ", fst_raw_max))
info(paste0("Windows with Fst > 1: ", fst_gt1_n))
info("Depth_ratio display max: 5")

# Circos
fst_ylim_circos   <- c(0, 1)
snps_ylim_circos  <- c(0, shared_snps_max)
depth_ylim_circos <- c(0, 5)

# -----------------------------
# Colors
# -----------------------------
COL_GREEN_DARK   <- "#7FC89A"
COL_GREEN_LIGHT  <- "#BFE6C9"

COL_PINK_DARK    <- "#E88AA3"
COL_PINK_LIGHT   <- "#F7B7C6"

COL_BLUE_DARK    <- "#7FB5E6"
COL_BLUE_LIGHT   <- "#B7D7F5"

COL_PURPLE_DARK  <- "#9D4EDD"
COL_PURPLE_LIGHT <- "#CDB4DB"

prefix <- detect_prefix()

# -----------------------------
# Manhattan
# -----------------------------
build_manhattan_plot(
  df_filtered = df_filt,
  keep_contigs = keep_contigs,
  shared_snps_max = shared_snps_max,
  fst_raw_max = fst_raw_max,
  fst_gt1_n = fst_gt1_n,
  prefix = prefix
)

# -----------------------------
# Circos tracks
# sem labels para evitar sobreposição superior do sgtr
# -----------------------------
tracks_circos <- list(
  single_metric_track(
    "Fst_circos",
    type = "points",
    label = "",
    colors = c(COL_GREEN_DARK, COL_GREEN_LIGHT),
    ylim = fst_ylim_circos
  ),
  single_metric_track(
    "Snps_females",
    type = "points",
    label = "",
    colors = c(COL_PINK_DARK, COL_PINK_LIGHT),
    ylim = snps_ylim_circos
  ),
  single_metric_track(
    "Snps_males",
    type = "points",
    label = "",
    colors = c(COL_BLUE_DARK, COL_BLUE_LIGHT),
    ylim = snps_ylim_circos
  ),
  single_metric_track(
    "Depth_ratio_circos",
    type = "points",
    label = "",
    colors = c(COL_PURPLE_DARK, COL_PURPLE_LIGHT),
    ylim = depth_ylim_circos
  )
)

info("Generating Circos plot with sgtr ...")

circos_tmp_png <- paste0(prefix, "_circos_FST_SNPf_SNPm_pastelAlt.tmp.png")
circos_out_png <- paste0(prefix, "_circos_FST_SNPf_SNPm_pastelAlt.png")

# PNG com legenda central manual
plot_circos(
  FILTERED_WIN,
  tracks = tracks_circos,
  chromosomes_file = CHR_TSV,
  output_file = circos_tmp_png
)

add_circos_center_legend_png(circos_tmp_png, circos_out_png)

if (file.exists(circos_tmp_png)) {
  unlink(circos_tmp_png)
}

info("DONE.")
info("Outputs:")
info(paste0("  - ", prefix, "_manhattan_FST_SNPf_SNPm_pastelAlt.png"))
info(paste0("  - ", prefix, "_circos_FST_SNPf_SNPm_pastelAlt.png"))
