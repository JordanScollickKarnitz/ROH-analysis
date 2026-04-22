# =============================================================================
# roh_metrics.R  —  Summary statistics derived from ROH segment tables
#
# FUNCTIONS:
#   compute_f_roh()           Genomic inbreeding coefficient per individual
#   summarize_roh_by_length() Partition ROH into temporal length classes
#   roh_summary_table()       Combined per-individual table for reporting
# =============================================================================


# =============================================================================
# compute_f_roh()
#
# The GENOMIC INBREEDING COEFFICIENT (F_ROH) is the fraction of the autosome
# inherited identical-by-descent from a recent common ancestor:
#
#   F_ROH = (total length of all ROH in individual) / (autosomal genome length)
#
# Interpretation:
#   F_ROH = 0.00  no detectable recent inbreeding
#   F_ROH = 0.25  equivalent to one parent being a full sibling of the other
#   F_ROH = 0.50  product of full-sibling x full-sibling mating
#
# Known approximate autosomal genome lengths for reference:
#   Elk (Cervus canadensis)    ~2.5 x 10^6 kb  (2.5 Gb)
#   White-tailed deer          ~2.7 x 10^6 kb
#   Gray wolf (Canis lupus)    ~2.4 x 10^6 kb
#   Human (Homo sapiens)       ~2.9 x 10^6 kb
#   Cattle (Bos taurus)        ~2.5 x 10^6 kb
#
# Args:
#   roh              data.frame of ROH segments (from detect_roh_*)
#   genome_length_kb total autosomal length in kilobases
#   all_individuals  optional character vector of ALL individual IDs —
#                    individuals with zero ROH will get F_ROH = 0 rather
#                    than being absent from the output
#
# Returns: data.frame with columns:
#   iid, n_roh, total_roh_kb, mean_roh_kb, max_roh_kb, f_roh
# =============================================================================
compute_f_roh <- function(roh, genome_length_kb, all_individuals = NULL) {

  if (missing(genome_length_kb) || is.null(genome_length_kb)) {
    stop(
      "You must supply genome_length_kb.\n",
      "For elk (Cervus canadensis), use: genome_length_kb = 2.5e6\n",
      "  (approximately 2.5 gigabases of autosomal sequence)"
    )
  }

  # Helper: compute stats for one individual's ROH
  .ind_stats <- function(d) {
    data.frame(
      iid          = unique(d$iid),
      n_roh        = nrow(d),
      total_roh_kb = sum(d$length_kb),
      mean_roh_kb  = mean(d$length_kb),
      max_roh_kb   = max(d$length_kb),
      f_roh        = sum(d$length_kb) / genome_length_kb,
      stringsAsFactors = FALSE
    )
  }

  if (nrow(roh) == 0) {
    base <- data.frame(iid = character(), n_roh = integer(),
                       total_roh_kb = numeric(), mean_roh_kb = numeric(),
                       max_roh_kb = numeric(), f_roh = numeric(),
                       stringsAsFactors = FALSE)
  } else {
    base <- do.call(rbind, lapply(split(roh, roh$iid), .ind_stats))
    rownames(base) <- NULL
  }

  # Fill in individuals with no ROH at all
  if (!is.null(all_individuals)) {
    missing_ids <- setdiff(all_individuals, base$iid)
    if (length(missing_ids) > 0) {
      zeros <- data.frame(
        iid          = missing_ids,
        n_roh        = 0L,
        total_roh_kb = 0,
        mean_roh_kb  = 0,
        max_roh_kb   = 0,
        f_roh        = 0,
        stringsAsFactors = FALSE
      )
      base <- rbind(base, zeros)
    }
  }

  base[order(-base$f_roh), ]
}


# =============================================================================
# summarize_roh_by_length()
#
# ROH length reflects how RECENT the inbreeding event was.
# Longer ROH = more recent common ancestor (less time for recombination to
# break up the shared haplotype).
#
# Standard length classes (based on Ceballos et al. 2018):
#
#   Class    Length (Mb)   Interpretation
#   ------   -----------   -----------------------------------------------
#   short     0.1 – 1      Very old relatedness (100+ generations)
#   medium    1   – 5      Moderately old (~20-50 generations)
#   long      5   – 16     Recent inbreeding (~5-25 generations)
#   vlong     > 16         Very recent close-kin mating (<5 generations)
#
# Args:
#   roh       data.frame of ROH segments
#   breaks_kb breakpoints in kilobases (default: 100, 1000, 5000, 16000, Inf)
#   labels    names for each class
#
# Returns: data.frame(iid, length_class, n, total_kb) — one row per
#          (individual x length class) combination that has at least one ROH.
# =============================================================================
summarize_roh_by_length <- function(roh,
                                     breaks_kb = c(0, 1000, 5000, 16000, Inf),
                                     labels    = c("short", "medium", "long", "vlong")) {

  if (nrow(roh) == 0) {
    return(data.frame(iid = character(), length_class = character(),
                      n = integer(), total_kb = numeric(),
                      stringsAsFactors = FALSE))
  }

  roh$length_class <- cut(roh$length_kb,
                           breaks = breaks_kb,
                           labels = labels,
                           include.lowest = TRUE)

  # Aggregate: count and sum per individual x class
  n_agg   <- stats::aggregate(length_kb ~ iid + length_class, data = roh, FUN = length)
  sum_agg <- stats::aggregate(length_kb ~ iid + length_class, data = roh, FUN = sum)

  colnames(n_agg)[3]   <- "n"
  colnames(sum_agg)[3] <- "total_kb"

  out <- merge(n_agg, sum_agg, by = c("iid", "length_class"))
  out[order(out$iid, out$length_class), ]
}


# =============================================================================
# roh_summary_table()
#
# Combines F_ROH and length-class breakdown into a single reporting table.
# Optionally merges in sample metadata (population, region, etc.) from fam.
# =============================================================================
roh_summary_table <- function(roh, genome_length_kb, fam = NULL) {

  base     <- compute_f_roh(roh, genome_length_kb = genome_length_kb,
                             all_individuals = if (!is.null(fam)) fam$iid else NULL)
  by_class <- summarize_roh_by_length(roh)

  # Pivot by_class to wide format (one column per length class)
  if (nrow(by_class) > 0) {
    wide <- stats::reshape(
      by_class[, c("iid", "length_class", "total_kb")],
      idvar     = "iid",
      timevar   = "length_class",
      direction = "wide"
    )
    # Rename: "total_kb.short" -> "kb_short"
    colnames(wide) <- gsub("^total_kb\\.", "kb_", colnames(wide))
    wide[is.na(wide)] <- 0
    out <- merge(base, wide, by = "iid", all.x = TRUE)
  } else {
    out <- base
  }

  # Attach population / other metadata
  if (!is.null(fam) && "iid" %in% colnames(fam)) {
    meta_cols <- setdiff(colnames(fam), colnames(out))
    meta_cols <- c("iid", meta_cols)
    out <- merge(fam[, intersect(meta_cols, colnames(fam))],
                 out, by = "iid", all.y = TRUE)
  }

  out[order(-out$f_roh), ]
}
