# =============================================================================
# plotting.R  —  ggplot2 visualizations of ROH analysis output
#
# All functions require ggplot2. Install with: install.packages("ggplot2")
# ggridges is optional (used in plot_roh_distribution for ridge plots).
#
# FUNCTIONS:
#   plot_f_roh_boxplot()     Boxplot of F_ROH by population
#   plot_roh_distribution()  Histogram / ridge plot of ROH lengths
#   plot_length_class_stack() Stacked bar of ROH kb by length class
#   plot_roh_karyotype()     Chromosome map of ROH segments per individual
#   plot_roh_island_density() Fraction of individuals with ROH per genomic bin
# =============================================================================

.check_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.\n",
         "Install with:  install.packages('ggplot2')")
  }
  invisible(TRUE)
}


# =============================================================================
# plot_f_roh_boxplot()
#
# Shows the distribution of individual genomic inbreeding (F_ROH) within
# each population. Points are jittered over the box to show raw data.
# Expected population order for elk tutorial: RANCH > TULE > ROCKY
# =============================================================================
plot_f_roh_boxplot <- function(f_roh_df, fam,
                                pop_order = NULL,
                                colors    = NULL,
                                title     = "Genomic inbreeding (F\u2080\u2080\u2080) by population") {
  .check_ggplot()

  d <- merge(f_roh_df[, c("iid", "f_roh")],
             fam[, c("iid", "population")], by = "iid")
  d <- d[!is.na(d$population), ]

  if (!is.null(pop_order)) {
    d$population <- factor(d$population, levels = pop_order)
  }

  default_colors <- c("ROCKY" = "#2166AC", "TULE" = "#F4A582", "RANCH" = "#D6604D")
  if (is.null(colors)) colors <- default_colors

  ggplot2::ggplot(d, ggplot2::aes(x = population, y = f_roh, fill = population)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    ggplot2::geom_jitter(ggplot2::aes(color = population),
                         width = 0.15, size = 1.8, alpha = 0.75) +
    ggplot2::scale_fill_manual(values  = colors, guide = "none") +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::labs(
      x     = NULL,
      y     = expression(F[ROH]),
      title = "Genomic inbreeding coefficient by elk population",
      subtitle = "Each point = one individual; boxes = median + IQR"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


# =============================================================================
# plot_roh_distribution()
#
# Histogram (or ridge plot) of ROH segment lengths, colored by population.
# X-axis is log10-scaled because ROH lengths span several orders of magnitude.
# =============================================================================
plot_roh_distribution <- function(roh, fam = NULL, by = "population",
                                   use_ridges = FALSE,
                                   title = "Distribution of ROH segment lengths") {
  .check_ggplot()

  if (!is.null(by) && !is.null(fam)) {
    if (by %in% colnames(fam)) {
      roh <- merge(roh, fam[, c("iid", by)], by = "iid", all.x = TRUE)
    }
  }

  roh$length_mb <- roh$length_kb / 1000

  if (use_ridges && requireNamespace("ggridges", quietly = TRUE) && !is.null(by)) {

    p <- ggplot2::ggplot(roh,
           ggplot2::aes(x = length_mb,
                        y = .data[[by]],
                        fill = .data[[by]])) +
      ggridges::geom_density_ridges(alpha = 0.65, scale = 1.2) +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::labs(x = "ROH length (Mb, log scale)", y = NULL,
                    title = title,
                    subtitle = "Longer ROH = more recent inbreeding")

  } else if (!is.null(by) && by %in% colnames(roh)) {

    p <- ggplot2::ggplot(roh,
           ggplot2::aes(x = length_mb, fill = .data[[by]])) +
      ggplot2::geom_histogram(alpha = 0.6, position = "identity",
                               bins = 40, color = NA) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "ROH length (Mb, log scale)", y = "Count",
                    fill = by, title = title)

  } else {

    p <- ggplot2::ggplot(roh, ggplot2::aes(x = length_mb)) +
      ggplot2::geom_histogram(bins = 40, fill = "steelblue4", color = NA) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "ROH length (Mb, log scale)", y = "Count", title = title)
  }

  p + ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


# =============================================================================
# plot_length_class_stack()
#
# Stacked bar chart showing total ROH (in Mb) per individual, partitioned
# into short / medium / long / vlong length classes.
# Individuals are grouped by population (facets).
#
# KEY INSIGHT: RANCH individuals should have tall bars dominated by "long" and
# "vlong" (blue/purple) segments. TULE individuals should have shorter bars
# dominated by "medium" (orange). ROCKY should have minimal bars.
# =============================================================================
plot_length_class_stack <- function(roh_length_summary, fam = NULL,
                                     title = "ROH length class composition per individual") {
  .check_ggplot()

  d <- roh_length_summary
  if (!is.null(fam)) {
    d <- merge(d, fam[, c("iid", "population")], by = "iid", all.x = TRUE)
  }
  d$total_mb <- d$total_kb / 1000

  # Color by biological meaning: red = very recent, blue = ancient
  class_colors <- c(
    "short"  = "#4393C3",   # blue — ancient
    "medium" = "#F4A582",   # orange — moderate
    "long"   = "#D6604D",   # red — recent
    "vlong"  = "#840006"    # dark red — very recent
  )

  p <- ggplot2::ggplot(d,
         ggplot2::aes(x = iid, y = total_mb, fill = length_class)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = class_colors,
                               name = "ROH class",
                               labels = c("short (<1 Mb, ancient)",
                                          "medium (1-5 Mb)",
                                          "long (5-16 Mb)",
                                          "vlong (>16 Mb, very recent)")) +
    ggplot2::labs(x = NULL, y = "Total ROH (Mb)", title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
      plot.title   = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(fam) && "population" %in% colnames(d)) {
    p <- p + ggplot2::facet_wrap(~ population, scales = "free_x", nrow = 1)
  }
  p
}


# =============================================================================
# plot_roh_karyotype()
#
# Chromosome-level view of ROH positions for a subset of individuals.
# Each horizontal bar is a ROH segment, colored by length.
# Useful for showing students what ROH actually look like on chromosomes.
# =============================================================================
plot_roh_karyotype <- function(roh, iids = NULL,
                                title = "ROH segments across chromosomes") {
  .check_ggplot()

  if (!is.null(iids)) roh <- roh[roh$iid %in% iids, ]
  if (nrow(roh) == 0) stop("No ROH to plot for the requested individuals.")

  roh$start_mb <- roh$start_bp / 1e6
  roh$end_mb   <- roh$end_bp   / 1e6

  ggplot2::ggplot(roh) +
    ggplot2::geom_segment(
      ggplot2::aes(x = start_mb, xend = end_mb,
                   y = iid, yend = iid,
                   color = length_kb / 1000),
      linewidth = 3.5, lineend = "round"
    ) +
    ggplot2::facet_wrap(~ paste("Chr", chr), scales = "free_x", nrow = 2) +
    ggplot2::scale_color_viridis_c(name = "ROH length\n(Mb)",
                                   option = "magma", direction = -1) +
    ggplot2::labs(x = "Position (Mb)", y = NULL, title = title,
                  subtitle = "Warmer colors = longer (more recent) ROH") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 7),
      plot.title   = ggplot2::element_text(face = "bold"),
      strip.text   = ggplot2::element_text(face = "bold")
    )
}


# =============================================================================
# plot_roh_island_density()
#
# Shows the fraction of individuals with a ROH at each position along one
# chromosome. Peaks = candidate ROH islands.
# =============================================================================
plot_roh_island_density <- function(roh, chromosome, bin_kb = 250,
                                     n_individuals = NULL,
                                     highlight_color = "#D6604D") {
  .check_ggplot()

  sub <- roh[roh$chr == as.character(chromosome), ]
  if (nrow(sub) == 0) stop("No ROH on chromosome: ", chromosome)
  if (is.null(n_individuals)) n_individuals <- length(unique(roh$iid))

  bin_bp  <- bin_kb * 1000L
  max_bp  <- max(sub$end_bp)
  n_bins  <- ceiling(max_bp / bin_bp)
  bin_mid <- ((seq_len(n_bins) - 0.5) * bin_bp) / 1e6

  frac <- vapply(seq_len(n_bins), function(b) {
    s <- (b - 1L) * bin_bp + 1L
    e <- b * bin_bp
    overlaps <- sub$start_bp <= e & sub$end_bp >= s
    length(unique(sub$iid[overlaps])) / n_individuals
  }, numeric(1))

  df <- data.frame(pos_mb = bin_mid, frac_ind = frac)

  thr <- stats::quantile(frac, 0.95)

  ggplot2::ggplot(df, ggplot2::aes(x = pos_mb, y = frac_ind)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = frac_ind >= thr),
      width = bin_kb / 1000 * 0.9
    ) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = highlight_color),
                               guide = "none") +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed", color = "grey30") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      x        = "Position (Mb)",
      y        = "Fraction of individuals with ROH",
      title    = paste0("ROH density — Chromosome ", chromosome),
      subtitle = paste0("Red bins = top 5% (candidate ROH islands); dashed = 95th percentile")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
