# =============================================================================
# pop_compare.R  —  Between-population comparisons of ROH burden
#
# FUNCTIONS:
#   compare_f_roh()          Test whether F_ROH differs across populations
#   compare_length_classes() Test whether ROH length composition differs
#   roh_islands()            Find genomic regions enriched for ROH
# =============================================================================


# =============================================================================
# compare_f_roh()
#
# Statistically compare genomic inbreeding (F_ROH) across populations.
#
# For two populations: Wilcoxon rank-sum test (non-parametric; F_ROH is
# typically not normally distributed).
# For three or more: Kruskal-Wallis test, followed by pairwise Wilcoxon
# post-hoc tests with Bonferroni correction.
#
# Args:
#   f_roh_df  data.frame with columns iid and f_roh (from compute_f_roh)
#   fam       data.frame with columns iid and population
#
# Returns: list with
#   $per_pop   summary stats per population (n, mean, median, sd)
#   $test      main test result (wilcox.test or kruskal.test object)
#   $posthoc   pairwise post-hoc results (only if > 2 populations)
# =============================================================================
compare_f_roh <- function(f_roh_df, fam) {

  stopifnot("iid"   %in% colnames(f_roh_df),
            "f_roh" %in% colnames(f_roh_df))
  stopifnot("iid"        %in% colnames(fam),
            "population" %in% colnames(fam))

  d <- merge(f_roh_df[, c("iid", "f_roh")],
             fam[, c("iid", "population")],
             by = "iid")
  d <- d[!is.na(d$population), ]

  # --- Per-population summary ---
  per_pop <- do.call(rbind, lapply(split(d, d$population), function(x) {
    data.frame(
      population   = unique(x$population),
      n            = nrow(x),
      mean_f_roh   = round(mean(x$f_roh,             na.rm = TRUE), 5),
      median_f_roh = round(stats::median(x$f_roh,    na.rm = TRUE), 5),
      sd_f_roh     = round(stats::sd(x$f_roh,        na.rm = TRUE), 5),
      min_f_roh    = round(min(x$f_roh,              na.rm = TRUE), 5),
      max_f_roh    = round(max(x$f_roh,              na.rm = TRUE), 5),
      stringsAsFactors = FALSE
    )
  }))
  rownames(per_pop) <- NULL
  per_pop <- per_pop[order(-per_pop$mean_f_roh), ]

  pops <- unique(d$population)

  # --- Main test ---
  if (length(pops) == 2) {
    grp1 <- d$f_roh[d$population == pops[1]]
    grp2 <- d$f_roh[d$population == pops[2]]
    main_test <- stats::wilcox.test(grp1, grp2,
                                     alternative = "two.sided",
                                     exact = FALSE)
    posthoc <- NULL
  } else if (length(pops) > 2) {
    main_test <- stats::kruskal.test(f_roh ~ as.factor(population), data = d)
    # Pairwise post-hoc with Bonferroni correction
    posthoc <- stats::pairwise.wilcox.test(
      x     = d$f_roh,
      g     = as.factor(d$population),
      p.adjust.method = "bonferroni",
      exact = FALSE
    )
  } else {
    main_test <- NULL
    posthoc   <- NULL
    message("  Only one population found — no test performed.")
  }

  list(per_pop = per_pop, test = main_test, posthoc = posthoc)
}


# =============================================================================
# compare_length_classes()
#
# Test whether populations differ in their ROH LENGTH DISTRIBUTION.
#
# A population with recent inbreeding has more LONG ROH.
# A population with an old bottleneck has more MEDIUM ROH.
# Comparing length classes tells you not just "who is more inbred" but
# "when did the inbreeding happen."
#
# Runs a Kruskal-Wallis test per length class on total kb per individual.
#
# Args:
#   roh_length_summary  output of summarize_roh_by_length()
#   fam                 data.frame with iid and population
#
# Returns: data.frame with test results per length class
# =============================================================================
compare_length_classes <- function(roh_length_summary, fam) {

  d <- merge(roh_length_summary, fam[, c("iid", "population")], by = "iid")
  d <- d[!is.na(d$population), ]

  classes <- levels(factor(d$length_class))

  out <- lapply(classes, function(cl) {
    sub <- d[d$length_class == cl, ]
    pops_present <- unique(sub$population)
    if (length(pops_present) < 2) return(NULL)

    kt <- stats::kruskal.test(total_kb ~ as.factor(population), data = sub)
    data.frame(
      length_class  = cl,
      H_statistic   = round(unname(kt$statistic), 3),
      df            = unname(kt$parameter),
      p_value       = round(kt$p.value, 4),
      interpretation = ifelse(kt$p.value < 0.05,
                              "Populations differ in this ROH class",
                              "No significant difference"),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), out))
}


# =============================================================================
# roh_islands()
#
# Identify "ROH islands" — genomic regions where an unusually high fraction
# of individuals carry a ROH. These can indicate:
#   - Regions under positive selection (selective sweep)
#   - Centromeres / low-recombination regions
#   - Recently shared common ancestors (population bottlenecks)
#
# Method: divide the genome into bins of `bin_kb` kilobases, count the
# fraction of individuals with a ROH overlapping each bin, flag bins above
# the `top_frac` quantile, merge adjacent flagged bins.
#
# Args:
#   roh           ROH segment data.frame
#   bin_kb        bin size in kilobases (default 500 = 0.5 Mb)
#   top_frac      quantile threshold for "elevated" (default 0.95)
#   n_individuals total number of individuals (denominator for frequency)
#
# Returns: data.frame of candidate ROH island regions, sorted by peak frequency
# =============================================================================
roh_islands <- function(roh, bin_kb = 500, top_frac = 0.95,
                         n_individuals = NULL) {

  if (nrow(roh) == 0) return(data.frame())
  if (is.null(n_individuals)) {
    n_individuals <- length(unique(roh$iid))
    message("  n_individuals not supplied — using n = ", n_individuals,
            " (individuals present in ROH table)")
  }

  bin_bp <- bin_kb * 1000L
  chroms <- unique(roh$chr)
  out    <- list()

  for (ch in chroms) {

    sub <- roh[roh$chr == ch, ]
    if (nrow(sub) == 0) next

    max_bp <- max(sub$end_bp)
    n_bins <- ceiling(max_bp / bin_bp)
    if (n_bins < 2) next

    bin_starts <- (seq_len(n_bins) - 1L) * bin_bp + 1L
    bin_ends   <- bin_starts + bin_bp - 1L

    # Count distinct individuals with a ROH overlapping each bin
    counts <- vapply(seq_len(n_bins), function(b) {
      overlaps <- sub$start_bp <= bin_ends[b] & sub$end_bp >= bin_starts[b]
      length(unique(sub$iid[overlaps]))
    }, integer(1))

    frac <- counts / n_individuals
    thr  <- stats::quantile(frac, top_frac, na.rm = TRUE)
    flag <- frac >= thr & frac > 0

    # Merge contiguous flagged bins
    runs  <- rle(flag)
    r_end <- cumsum(runs$lengths)
    r_str <- r_end - runs$lengths + 1L

    for (k in which(runs$values)) {
      s <- r_str[k]; e <- r_end[k]
      out[[length(out) + 1L]] <- data.frame(
        chr           = ch,
        start_bp      = bin_starts[s],
        end_bp        = bin_ends[e],
        n_bins        = e - s + 1L,
        length_mb     = round((bin_ends[e] - bin_starts[s]) / 1e6, 2),
        peak_frac     = round(max(frac[s:e]), 3),
        mean_frac     = round(mean(frac[s:e]), 3),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(data.frame())
  res <- do.call(rbind, out)
  res[order(-res$peak_frac, res$chr), ]
}
