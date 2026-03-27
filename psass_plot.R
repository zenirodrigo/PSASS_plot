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
      "Please run this script inside the PSASS analyze folder."
    ))
  }
}

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

parse_chr_from_cli <- function() {
  v <- get_arg_value("--chr")
  if (is.na(v)) return(NA_integer_)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 1) {
    die("Invalid value for --chr. Use a positive integer, e.g. --chr 11")
  }
  n
}

parse_region_from_cli <- function() {
  v <- get_arg_value("--region")
  if (is.na(v)) return(NULL)

  if (!grepl("^[0-9]+:[0-9]+$", v)) {
    die("Invalid value for --region. Use start:end, e.g. --region 1:100000")
  }

  parts <- strsplit(v, ":", fixed = TRUE)[[1]]
  start <- suppressWarnings(as.numeric(parts[1]))
  end   <- suppressWarnings(as.numeric(parts[2]))

  if (is.na(start) || is.na(end) || start < 1 || end < 1 || end < start) {
    die("Invalid region coordinates. Use positive integers with end >= start, e.g. --region 1:100000")
  }

  list(start = start, end = end, raw = v)
}

parse_p_from_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  "--p" %in% args
}

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

detect_prefix <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (grepl("psass_hap1", wd, ignore.case = TRUE)) return("HAP1_PSASS")
  if (grepl("psass_hap2", wd, ignore.case = TRUE)) return("HAP2_PSASS")
  return("PSASS")
}

order_contigs_obvious <- function(contig_sizes_df) {
  contigs <- contig_sizes_df$Contig

  extract_chr_number <- function(x) {
    m <- regexec("^(Chromosome|chromosome|Chr|chr)[_\\-]?([0-9]+)$", x)
    reg <- regmatches(x, m)
    out <- rep(NA_integer_, length(x))

    for (i in seq_along(reg)) {
      if (length(reg[[i]]) >= 3) {
        out[i] <- suppressWarnings(as.integer(reg[[i]][3]))
      }
    }
    out
  }

  chr_num <- extract_chr_number(contigs)

  if (any(!is.na(chr_num))) {
    numeric_contigs <- contigs[!is.na(chr_num)]
    numeric_nums    <- chr_num[!is.na(chr_num)]
    other_contigs   <- contigs[is.na(chr_num)]

    numeric_ordered <- numeric_contigs[order(numeric_nums)]
    other_ordered   <- sort(other_contigs)

    return(c(numeric_ordered, other_ordered))
  }

  extract_scf_number <- function(x) {
    m <- regexec("^(scaffold|Scaffold)[_\\-]?([0-9]+)$", x)
    reg <- regmatches(x, m)
    out <- rep(NA_integer_, length(x))

    for (i in seq_along(reg)) {
      if (length(reg[[i]]) >= 3) {
        out[i] <- suppressWarnings(as.integer(reg[[i]][3]))
      }
    }
    out
  }

  sca_num <- extract_scf_number(contigs)

  if (any(!is.na(sca_num))) {
    numeric_contigs <- contigs[!is.na(sca_num)]
    numeric_nums    <- sca_num[!is.na(sca_num)]
    other_contigs   <- contigs[is.na(sca_num)]

    numeric_ordered <- numeric_contigs[order(numeric_nums)]
    other_ordered   <- sort(other_contigs)

    return(c(numeric_ordered, other_ordered))
  }

  contig_sizes_df <- contig_sizes_df[order(contig_sizes_df$Length, decreasing = TRUE), ]
  contig_sizes_df$Contig
}

write_chromosomes_tsv <- function(contigs_ordered, outfile) {
  labels <- contigs_ordered
  labels <- sub("^Chromosome", "Chr", labels)
  labels <- sub("^(scaffold|Scaffold)[_\\-]?", "Scf", labels)

  out <- data.frame(Contig = contigs_ordered, Label = labels, stringsAsFactors = FALSE)
  write.table(out, outfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  info(paste0("Wrote chromosomes mapping file: ", outfile))
}

sanitize_for_filename <- function(x) {
  x <- gsub(":", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x
}

# -----------------------------
# X-axis helper for --chr / --chr --region
# -----------------------------
build_single_chr_x_axis <- function(region_info = NULL, chromosome_length_bp = NULL) {

  if (!is.null(region_info)) {
    span_bp <- region_info$end - region_info$start

    if (span_bp < 1e6) {
      unit_bp <- 1e3
      unit_label <- "kb"
    } else {
      unit_bp <- 1e6
      unit_label <- "Mb"
    }

    major_step_bp <- 5 * unit_bp
    minor_step_bp <- 1 * unit_bp

    axis_start_bp <- floor(region_info$start / major_step_bp) * major_step_bp
    axis_end_bp   <- ceiling(region_info$end / major_step_bp) * major_step_bp

  } else {
    if (is.null(chromosome_length_bp) || !is.finite(chromosome_length_bp) || chromosome_length_bp <= 0) {
      die("Could not determine chromosome length for single-chromosome Manhattan axis.")
    }

    unit_bp <- 1e6
    unit_label <- "Mb"
    major_step_bp <- 5 * unit_bp
    minor_step_bp <- 1 * unit_bp

    axis_start_bp <- 0
    axis_end_bp   <- ceiling(chromosome_length_bp / major_step_bp) * major_step_bp
  }

  breaks_bp <- seq(axis_start_bp, axis_end_bp, by = minor_step_bp)
  labels <- ifelse((breaks_bp %% major_step_bp) == 0,
                   format(round(breaks_bp / unit_bp), trim = TRUE, scientific = FALSE),
                   "")

  list(
    unit_bp = unit_bp,
    unit_label = unit_label,
    axis_start_bp = axis_start_bp,
    axis_end_bp = axis_end_bp,
    breaks_scaled = breaks_bp / unit_bp,
    labels = labels,
    limits_scaled = c(axis_start_bp, axis_end_bp) / unit_bp
  )
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
# Convert PNG to PDF
# -----------------------------
png_to_pdf <- function(input_png, output_pdf) {
  img <- png::readPNG(input_png)
  h <- dim(img)[1]
  w <- dim(img)[2]

  grDevices::pdf(output_pdf, width = w / 150, height = h / 150, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.raster(img, x = 0.5, y = 0.5, width = 1, height = 1)
  grDevices::dev.off()

  info(paste0("Wrote PDF version: ", output_pdf))
}

# -----------------------------
# Manhattan helpers
# -----------------------------
build_manhattan_plot <- function(df_filtered,
                                 keep_contigs,
                                 shared_snps_max,
                                 fst_raw_max,
                                 fst_gt1_n,
                                 prefix,
                                 subtitle_extra = NULL,
                                 single_chr_axis = FALSE,
                                 region_info = NULL,
                                 fill_under_points = FALSE) {

  df_plot <- df_filtered %>%
    mutate(
      Position_Mbp = Position / 1e6,
      Depth_ratio_plot = pmin(Depth_ratio, 5)
    )

  chrom_lengths <- df_plot %>%
    group_by(Chromosome = Contig) %>%
    summarise(
      Length_bp = max(Length, na.rm = TRUE),
      Length_Mbp = max(Length, na.rm = TRUE) / 1e6,
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
    ) %>%
    arrange(factor(Contig, levels = keep_contigs), Position)

  y_max_snps <- shared_snps_max * 1.05
  y_max_fst  <- min(fst_raw_max * 1.05, 1)
  y_max_depth_ratio <- 5

  chr_labels <- gsub("^Chromosome", "Chr", chrom_lengths$Chromosome)
  chr_labels <- gsub("^(scaffold|Scaffold)[_\\-]?", "Scf", chr_labels)

  info("Generating Manhattan plot with ggplot2 ...")
  info(paste0("  Manhattan SNP Y max (shared): ", round(y_max_snps, 2)))
  info(paste0("  Manhattan Fst display max: ", round(y_max_fst, 4)))
  info("  Manhattan Depth_ratio display max: 5")

  x_scale_single <- NULL
  if (single_chr_axis) {
    chr_len_bp <- max(df_plot$Length, na.rm = TRUE)

    x_scale_single <- build_single_chr_x_axis(
      region_info = region_info,
      chromosome_length_bp = chr_len_bp
    )

    df_plot <- df_plot %>%
      mutate(
        X_axis_scaled = Position / x_scale_single$unit_bp
      )

    info(paste0(
      "  Single-chromosome X axis enabled in ",
      x_scale_single$unit_label,
      " | start=",
      round(x_scale_single$limits_scaled[1], 3),
      " | end=",
      round(x_scale_single$limits_scaled[2], 3)
    ))
  }

  create_male_plot <- function(data, y_max) {
    p <- ggplot(data)

    if (single_chr_axis) {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = X_axis_scaled, y = Snps_males),
            fill = COL_BLUE_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
        geom_point(
          aes(x = X_axis_scaled, y = Snps_males),
          color = COL_BLUE_DARK,
          size = 0.3,
          alpha = 0.8,
          shape = 16
        ) +
        geom_smooth(
          aes(x = X_axis_scaled, y = Snps_males),
          method = "loess",
          se = FALSE,
          color = COL_BLUE_DARK,
          linewidth = 1.0,
          span = 0.05
        ) +
        scale_x_continuous(
          name = NULL,
          breaks = x_scale_single$breaks_scaled,
          labels = x_scale_single$labels,
          limits = x_scale_single$limits_scaled,
          expand = c(0, 0)
        )
    } else {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = Cumulative_position_Mbp, y = Snps_males, group = 1),
            fill = COL_BLUE_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
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
        )
    }

    p +
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
        axis.text.x = if (single_chr_axis) element_text(size = 8, color = "black") else element_blank(),
        axis.ticks.x = if (single_chr_axis) element_line(color = "black", linewidth = 0.3) else element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(5, 10, 0, 10)
      )
  }

  create_female_plot <- function(data, y_max) {
    p <- ggplot(data)

    if (single_chr_axis) {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = X_axis_scaled, y = Snps_females),
            fill = COL_PINK_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
        geom_point(
          aes(x = X_axis_scaled, y = Snps_females),
          color = COL_PINK_DARK,
          size = 0.3,
          alpha = 0.8,
          shape = 16
        ) +
        geom_smooth(
          aes(x = X_axis_scaled, y = Snps_females),
          method = "loess",
          se = FALSE,
          color = COL_PINK_DARK,
          linewidth = 1.0,
          span = 0.05
        ) +
        scale_x_continuous(
          name = NULL,
          breaks = x_scale_single$breaks_scaled,
          labels = x_scale_single$labels,
          limits = x_scale_single$limits_scaled,
          expand = c(0, 0)
        )
    } else {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = Cumulative_position_Mbp, y = Snps_females, group = 1),
            fill = COL_PINK_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
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
        )
    }

    p +
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
        axis.text.x = if (single_chr_axis) element_text(size = 8, color = "black") else element_blank(),
        axis.ticks.x = if (single_chr_axis) element_line(color = "black", linewidth = 0.3) else element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(0, 10, 0, 10)
      )
  }

  create_fst_plot <- function(data, y_max) {
    p <- ggplot(data)

    if (single_chr_axis) {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = X_axis_scaled, y = Fst),
            fill = COL_GREEN_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
        geom_point(
          aes(x = X_axis_scaled, y = Fst),
          color = COL_GREEN_DARK,
          size = 0.3,
          alpha = 0.8,
          shape = 16
        ) +
        geom_smooth(
          aes(x = X_axis_scaled, y = Fst),
          method = "loess",
          se = FALSE,
          color = COL_GREEN_DARK,
          linewidth = 1.0,
          span = 0.05
        ) +
        scale_x_continuous(
          name = NULL,
          breaks = x_scale_single$breaks_scaled,
          labels = x_scale_single$labels,
          limits = x_scale_single$limits_scaled,
          expand = c(0, 0)
        )
    } else {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = Cumulative_position_Mbp, y = Fst, group = 1),
            fill = COL_GREEN_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
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
        scale_color_manual(
          values = c("odd" = COL_GREEN_LIGHT, "even" = COL_GREEN_DARK),
          guide = "none"
        ) +
        scale_x_continuous(
          name = NULL,
          breaks = chrom_lengths$Midpoint_Mbp,
          labels = chr_labels,
          expand = expansion(mult = 0.02)
        )
    }

    p +
      geom_hline(
        yintercept = 0.25,
        linetype = "dashed",
        color = "gray30",
        alpha = 0.8
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
        axis.text.x = if (single_chr_axis) element_text(size = 8, color = "black") else element_blank(),
        axis.ticks.x = if (single_chr_axis) element_line(color = "black", linewidth = 0.3) else element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(0, 10, 0, 10)
      )
  }

  create_depth_ratio_plot <- function(data, y_max) {
    p <- ggplot(data)

    if (single_chr_axis) {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = X_axis_scaled, y = Depth_ratio_plot),
            fill = COL_PURPLE_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
        geom_point(
          aes(x = X_axis_scaled, y = Depth_ratio_plot),
          color = COL_PURPLE_DARK,
          size = 0.3,
          alpha = 0.8,
          shape = 16
        ) +
        geom_smooth(
          aes(x = X_axis_scaled, y = Depth_ratio_plot),
          method = "loess",
          se = FALSE,
          color = COL_PURPLE_DARK,
          linewidth = 1.0,
          span = 0.05
        ) +
        scale_x_continuous(
          name = paste0("Position (", x_scale_single$unit_label, ")"),
          breaks = x_scale_single$breaks_scaled,
          labels = x_scale_single$labels,
          limits = x_scale_single$limits_scaled,
          expand = c(0, 0)
        )
    } else {
      if (fill_under_points) {
        p <- p +
          geom_area(
            aes(x = Cumulative_position_Mbp, y = Depth_ratio_plot, group = 1),
            fill = COL_PURPLE_LIGHT,
            alpha = 0.35,
            linewidth = 0
          )
      }

      p <- p +
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
        scale_color_manual(
          values = c("odd" = COL_PURPLE_LIGHT, "even" = COL_PURPLE_DARK),
          guide = "none"
        ) +
        scale_x_continuous(
          name = "Chromosome",
          breaks = chrom_lengths$Midpoint_Mbp,
          labels = chr_labels,
          expand = expansion(mult = 0.02)
        )
    }

    p +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "gray30",
        alpha = 0.8
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
        axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, color = "black"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
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

  subtitle_text <- paste(
    "Scaffolds:", length(keep_contigs),
    "| Same Y scale for SNPs (max:", round(shared_snps_max, 0), ")"
  )

  if (!is.null(subtitle_extra) && nzchar(subtitle_extra)) {
    subtitle_text <- paste(subtitle_text, "|", subtitle_extra)
  }

  if (fill_under_points) {
    subtitle_text <- paste(subtitle_text, "|", "Mode: --p enabled")
  }

  combined_plot <- p_male / p_female / p_fst / p_depth +
    plot_layout(heights = c(1, 1, 1.2, 1.1)) +
    plot_annotation(
      title = "PSASS Analysis: SNP counts, Fst distribution and depth ratio",
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )

  out_png <- paste0(prefix, "_manhattan_FST_SNPf_SNPm_pastelAlt.png")
  out_pdf <- paste0(prefix, "_manhattan_FST_SNPf_SNPm_pastelAlt.pdf")

  ggsave(
    out_png,
    combined_plot,
    width = 16,
    height = 17,
    device = "png",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    out_pdf,
    combined_plot,
    width = 16,
    height = 17,
    device = cairo_pdf,
    bg = "white"
  )

  info(paste0("Wrote Manhattan plot PNG: ", out_png))
  info(paste0("Wrote Manhattan plot PDF: ", out_pdf))
}


REQ_WINDOW   <- "psass_window.tsv"
CLEAN_WINDOW <- "psass_window.clean.tsv"
FILTERED_WIN <- "psass_window.filtered.tsv"
CHR_TSV      <- "chromosomes.tsv"
SELECTED_TXT <- "selected_contigs.txt"

REQUESTED_N      <- parse_n_from_cli()
REQUESTED_CHR    <- parse_chr_from_cli()
REQUESTED_REGION <- parse_region_from_cli()
REQUESTED_P      <- parse_p_from_cli()

info("------------------------------------------------------------")
info("PSASS plotting helper")
info("Input required in this folder: psass_window.tsv")
info("Usage examples:")
info("  Rscript PSASS_plot.R --n 25")
info("  Rscript PSASS_plot.R --chr 11")
info("  Rscript PSASS_plot.R --chr 11 --region 1:100000")
info("  Rscript PSASS_plot.R --p")
info("  Rscript PSASS_plot.R --chr 11 --p")
info("  Rscript PSASS_plot.R --chr 11 --region 1:100000 --p")
info("  Rscript PSASS_plot.R")
info("------------------------------------------------------------")

if (!is.na(REQUESTED_N) && !is.na(REQUESTED_CHR)) {
  die("Use either --n OR --chr, not both together.")
}

if (!is.null(REQUESTED_REGION) && is.na(REQUESTED_CHR)) {
  die("The option --region requires --chr. Example: --chr 11 --region 1:100000")
}

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

df$Position      <- suppressWarnings(as.numeric(df$Position))
df$Fst           <- suppressWarnings(as.numeric(df$Fst))
df$Snps_females  <- suppressWarnings(as.numeric(df$Snps_females))
df$Snps_males    <- suppressWarnings(as.numeric(df$Snps_males))
df$Length        <- suppressWarnings(as.numeric(df$Length))
df$Depth_ratio   <- suppressWarnings(as.numeric(df$Depth_ratio))

contig_sizes <- aggregate(Length ~ Contig, data = df, FUN = max)
ordered_contigs <- order_contigs_obvious(contig_sizes)
max_n <- length(ordered_contigs)

region_tag <- NULL
subtitle_extra <- NULL

if (!is.na(REQUESTED_CHR)) {
  if (REQUESTED_CHR > max_n) {
    die(paste0(
      "Invalid value for --chr: ", REQUESTED_CHR,
      ". There are only ", max_n, " sequences/contigs available."
    ))
  }

  selected_contig <- ordered_contigs[REQUESTED_CHR]
  keep_contigs <- selected_contig

  info(paste0(
    "Using only sequence index ", REQUESTED_CHR,
    " -> contig: ", selected_contig
  ))

  if (!is.null(REQUESTED_REGION)) {
    region_tag <- sanitize_for_filename(paste0("region_", REQUESTED_REGION$raw))
    subtitle_extra <- paste0(
      "Contig index: ", REQUESTED_CHR,
      " (", selected_contig, ")",
      " | Region: ", REQUESTED_REGION$start, "-", REQUESTED_REGION$end
    )
  } else {
    subtitle_extra <- paste0(
      "Contig index: ", REQUESTED_CHR,
      " (", selected_contig, ")"
    )
  }

} else if (!is.na(REQUESTED_N)) {
  n_keep <- min(REQUESTED_N, max_n)
  if (REQUESTED_N > max_n) {
    info(paste0("You requested ", REQUESTED_N, " contigs, but only ", max_n, " are available. Using ALL (", max_n, ")."))
  } else {
    info(paste0("Keeping the first ", n_keep, " contigs for analysis."))
  }
  keep_contigs <- head(ordered_contigs, n_keep)
  subtitle_extra <- paste0("Top contigs used: ", length(keep_contigs))

} else {
  keep_contigs <- ordered_contigs
  info(paste0("No --n or --chr provided. Keeping ALL contigs for analysis (", max_n, ")."))
  info("Tip: if your assembly is fragmented (many scaffolds), use e.g. --n 25 or --n 50 to avoid Circos errors.")
  subtitle_extra <- paste0("All contigs used: ", length(keep_contigs))
}

writeLines(keep_contigs, SELECTED_TXT)
info(paste0("Wrote selected contigs list: ", SELECTED_TXT))

df_filt <- df[df$Contig %in% keep_contigs, ]

if (nrow(df_filt) == 0) {
  die("After filtering contigs, no rows remained in the dataset.")
}

if (!is.null(REQUESTED_REGION)) {
  chosen_contig <- keep_contigs[1]

  df_filt <- df_filt[
    df_filt$Contig == chosen_contig &
      df_filt$Position >= REQUESTED_REGION$start &
      df_filt$Position <= REQUESTED_REGION$end,
  ]

  if (nrow(df_filt) == 0) {
    die(paste0(
      "No rows remained after applying region filter ",
      REQUESTED_REGION$start, ":", REQUESTED_REGION$end,
      " on contig ", chosen_contig
    ))
  }

  info(paste0(
    "Applied region filter on ", chosen_contig,
    ": ", REQUESTED_REGION$start, "-", REQUESTED_REGION$end
  ))
}

df_filt$Fst_circos <- pmax(pmin(df_filt$Fst, 1), 0)
df_filt$Depth_ratio_circos <- pmin(df_filt$Depth_ratio, 5)

write.table(df_filt, FILTERED_WIN, sep = "\t", row.names = FALSE, quote = FALSE)
info(paste0("Wrote filtered window file: ", FILTERED_WIN))

write_chromosomes_tsv(keep_contigs, CHR_TSV)

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

fst_ylim_circos   <- c(0, 1)
snps_ylim_circos  <- c(0, shared_snps_max)
depth_ylim_circos <- c(0, 5)

COL_GREEN_DARK   <- "#7FC89A"
COL_GREEN_LIGHT  <- "#BFE6C9"

COL_PINK_DARK    <- "#E88AA3"
COL_PINK_LIGHT   <- "#F7B7C6"

COL_BLUE_DARK    <- "#7FB5E6"
COL_BLUE_LIGHT   <- "#B7D7F5"

COL_PURPLE_DARK  <- "#9D4EDD"
COL_PURPLE_LIGHT <- "#CDB4DB"

prefix_base <- detect_prefix()

if (!is.na(REQUESTED_CHR)) {
  prefix_base <- paste0(prefix_base, "_chr", REQUESTED_CHR)
}
if (!is.null(region_tag)) {
  prefix_base <- paste0(prefix_base, "_", region_tag)
}

build_manhattan_plot(
  df_filtered = df_filt,
  keep_contigs = keep_contigs,
  shared_snps_max = shared_snps_max,
  fst_raw_max = fst_raw_max,
  fst_gt1_n = fst_gt1_n,
  prefix = prefix_base,
  subtitle_extra = subtitle_extra,
  single_chr_axis = (!is.na(REQUESTED_CHR)),
  region_info = REQUESTED_REGION,
  fill_under_points = REQUESTED_P
)

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

circos_tmp_png <- paste0(prefix_base, "_circos_FST_SNPf_SNPm_pastelAlt.tmp.png")
circos_out_png <- paste0(prefix_base, "_circos_FST_SNPf_SNPm_pastelAlt.png")
circos_out_pdf <- paste0(prefix_base, "_circos_FST_SNPf_SNPm_pastelAlt.pdf")

plot_circos(
  FILTERED_WIN,
  tracks = tracks_circos,
  chromosomes_file = CHR_TSV,
  output_file = circos_tmp_png
)

add_circos_center_legend_png(circos_tmp_png, circos_out_png)
png_to_pdf(circos_out_png, circos_out_pdf)

if (file.exists(circos_tmp_png)) {
  unlink(circos_tmp_png)
}

info("DONE.")
info("Correct usage example:")
info("  Rscript PSASS_plot.R --chr 2 --region 50000:350000")
info("  Rscript PSASS_plot.R --p")
info("  Rscript PSASS_plot.R --chr 2 --region 50000:350000 --p")
info("Outputs:")
info(paste0("  - ", prefix_base, "_manhattan_FST_SNPf_SNPm_pastelAlt.png"))
info(paste0("  - ", prefix_base, "_manhattan_FST_SNPf_SNPm_pastelAlt.pdf"))
info(paste0("  - ", prefix_base, "_circos_FST_SNPf_SNPm_pastelAlt.png"))
info(paste0("  - ", prefix_base, "_circos_FST_SNPf_SNPm_pastelAlt.pdf"))
