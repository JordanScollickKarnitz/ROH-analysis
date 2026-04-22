# =============================================================================
# roh_detection.R  —  Algorithms for calling runs of homozygosity (ROH)
#
# A run of homozygosity is a long stretch of consecutive SNPs where an
# individual carries identical alleles on both chromosome copies (homozygous).
# These arise when both chromosomes were inherited from a common ancestor
# (identical-by-descent, IBD).
#
# TWO METHODS:
#
#   1. detect_roh_sliding()     — PLINK-style sliding-window algorithm
#      The standard in the field. Uses overlapping windows to tolerate
#      occasional genotyping errors (allows a few hets per window).
#      Best choice for real data with any genotyping error.
#
#   2. detect_roh_consecutive() — Consecutive homozygous SNP runs
#      Simpler and more transparent. Any single het or missing call breaks
#      the run (unless allow_missing = TRUE). Better for teaching the concept;
#      less robust to error in real data.
#
# INTERNAL GENOTYPE CODING (must be 0/1/2/NA throughout):
#   0  =  homozygous reference   (AA)
#   1  =  heterozygous           (Aa)
#   2  =  homozygous alternate   (aa)
#   NA =  missing / failed call
#
# INPUT: the list returned by any io_*.R reader:
#   geno$genotypes  integer matrix [n_ind x n_snp]
#   geno$map        data.frame(chr, snp_id, cm, bp)
#   geno$fam        data.frame(fid, iid, ..., population)
#
# OUTPUT: data.frame, one row per ROH segment, with columns:
#   iid, chr, start_bp, end_bp, length_kb, n_snps, n_het, n_missing
# =============================================================================


# ---- Internal helper -------------------------------------------------------
# Returns TRUE where a genotype is homozygous (0 or 2) and not NA
.is_hom <- function(g) !is.na(g) & (g == 0L | g == 2L)


# =============================================================================
# detect_roh_sliding()
#
# How the algorithm works (step by step):
#
#   For each individual, for each chromosome:
#
#   Step 1 — Slide a window of `window_snps` SNPs across the chromosome.
#             For each window position, check:
#               Is the number of HET genotypes <= max_het_in_window?
#               Is the number of MISSING genotypes <= max_missing_in_window?
#             If both conditions pass, label this a "homozygous window."
#
#   Step 2 — For each SNP, count how many of the windows overlapping it
#             were labeled "homozygous." Compute the fraction.
#
#   Step 3 — Label a SNP as "in ROH" if:
#               fraction >= window_threshold_frac   (default 5%)
#               AND the SNP itself is homozygous
#
#   Step 4 — Find consecutive stretches of "in ROH" SNPs.
#             Keep stretches that are >= min_snps long AND >= min_length_kb.
#
# Parameters:
#   window_snps           - size of the sliding window in SNPs (default 50)
#   min_snps              - minimum ROH length in SNPs (default 100)
#   min_length_kb         - minimum ROH length in kilobases (default 1000 = 1 Mb)
#   max_het_in_window     - heterozygous SNPs tolerated per window (default 1)
#   max_missing_in_window - missing SNPs tolerated per window (default 2)
#   window_threshold_frac - minimum fraction of overlapping windows that must
#                           be homozygous for a SNP to be called "in ROH"
# =============================================================================
detect_roh_sliding <- function(geno,
                                window_snps           = 50,
                                min_snps              = 100,
                                min_length_kb         = 1000,
                                max_het_in_window     = 1,
                                max_missing_in_window = 2,
                                window_threshold_frac = 0.05,
                                verbose               = TRUE) {

  stopifnot(is.list(geno), is.matrix(geno$genotypes))
  stopifnot(all(c("chr", "bp") %in% colnames(geno$map)))

  G    <- geno$genotypes
  map  <- geno$map
  iids <- if (!is.null(geno$fam$iid)) geno$fam$iid else rownames(G)
  if (is.null(iids)) iids <- paste0("IND", seq_len(nrow(G)))

  chroms   <- unique(map$chr)
  roh_list <- list()

  for (ind_i in seq_len(nrow(G))) {

    if (verbose && (ind_i == 1 || ind_i %% 10 == 0))
      message(sprintf("  Scanning individual %d / %d  [%s]",
                      ind_i, nrow(G), iids[ind_i]))

    for (ch in chroms) {

      # --- Pull out this chromosome's SNP indices and positions ---
      idx   <- which(map$chr == ch)
      if (length(idx) < window_snps) next

      gvec  <- G[ind_i, idx]           # genotypes for this individual on this chr
      bpvec <- map$bp[idx]             # base-pair positions

      # --- Per-SNP boolean flags ---
      hom_snp <- .is_hom(gvec)
      het_snp <- !is.na(gvec) & gvec == 1L
      mis_snp <- is.na(gvec)

      n_snp_chr <- length(gvec)
      n_windows  <- n_snp_chr - window_snps + 1L
      if (n_windows < 1L) next

      # --- Count how many overlapping windows each SNP falls in,
      #     and how many of those were "homozygous windows" ---
      snp_window_hits  <- integer(n_snp_chr)
      snp_window_total <- integer(n_snp_chr)

      # Rolling totals — initialize on the first window
      het_count <- sum(het_snp[seq_len(window_snps)])
      mis_count <- sum(mis_snp[seq_len(window_snps)])

      for (w in seq_len(n_windows)) {

        win_idx    <- w:(w + window_snps - 1L)
        is_hom_win <- (het_count <= max_het_in_window) &&
                      (mis_count <= max_missing_in_window)

        snp_window_total[win_idx] <- snp_window_total[win_idx] + 1L
        if (is_hom_win) {
          snp_window_hits[win_idx] <- snp_window_hits[win_idx] + 1L
        }

        # Slide the window forward by 1 SNP:
        #   subtract the SNP leaving on the left, add the SNP entering on the right
        if (w < n_windows) {
          het_count <- het_count - het_snp[w] + het_snp[w + window_snps]
          mis_count <- mis_count - mis_snp[w] + mis_snp[w + window_snps]
        }
      }

      # --- Flag SNPs as "in ROH" ---
      frac   <- ifelse(snp_window_total > 0,
                       snp_window_hits / snp_window_total,
                       0)
      in_roh <- frac >= window_threshold_frac & hom_snp

      # --- Find consecutive runs of in_roh == TRUE ---
      runs  <- rle(in_roh)
      r_end <- cumsum(runs$lengths)
      r_str <- r_end - runs$lengths + 1L

      for (k in which(runs$values & runs$lengths >= min_snps)) {
        s          <- r_str[k]
        e          <- r_end[k]
        length_kb  <- (bpvec[e] - bpvec[s]) / 1000

        if (length_kb < min_length_kb) next

        seg <- gvec[s:e]
        roh_list[[length(roh_list) + 1L]] <- data.frame(
          iid       = iids[ind_i],
          chr       = ch,
          start_bp  = bpvec[s],
          end_bp    = bpvec[e],
          length_kb = length_kb,
          n_snps    = e - s + 1L,
          n_het     = sum(!is.na(seg) & seg == 1L),
          n_missing = sum(is.na(seg)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(roh_list) == 0) {
    if (verbose) message("  No ROH segments passed filters.")
    return(.empty_roh_df())
  }

  result <- do.call(rbind, roh_list)
  if (verbose) message(sprintf("  Done: %d ROH segments found", nrow(result)))
  result
}


# =============================================================================
# detect_roh_consecutive()
#
# Simpler algorithm — every heterozygous or missing call breaks the run.
#
# Steps:
#   1. For each individual x chromosome, label each SNP TRUE (homozygous)
#      or FALSE (het / missing).
#   2. Find maximal consecutive TRUE stretches.
#   3. Keep stretches >= min_snps AND >= min_length_kb.
#
# Setting allow_missing = TRUE lets missing calls pass through without
# breaking the run (useful when missingness is very low and random).
#
# This method is more sensitive to genotyping error than detect_roh_sliding().
# Use it to understand the algorithm, then switch to the sliding-window
# method for your actual analysis.
# =============================================================================
detect_roh_consecutive <- function(geno,
                                    min_snps      = 50,
                                    min_length_kb = 500,
                                    allow_missing = FALSE,
                                    verbose       = TRUE) {

  G    <- geno$genotypes
  map  <- geno$map
  iids <- if (!is.null(geno$fam$iid)) geno$fam$iid else rownames(G)
  if (is.null(iids)) iids <- paste0("IND", seq_len(nrow(G)))

  chroms   <- unique(map$chr)
  roh_list <- list()

  for (ind_i in seq_len(nrow(G))) {

    if (verbose && (ind_i == 1 || ind_i %% 10 == 0))
      message(sprintf("  Scanning individual %d / %d  [%s]",
                      ind_i, nrow(G), iids[ind_i]))

    for (ch in chroms) {

      idx   <- which(map$chr == ch)
      if (length(idx) < min_snps) next

      gvec  <- G[ind_i, idx]
      bpvec <- map$bp[idx]

      # Label each SNP as "passable" (can contribute to a ROH)
      passable <- if (allow_missing) {
        .is_hom(gvec) | is.na(gvec)
      } else {
        .is_hom(gvec)   # het or missing both break the run
      }

      runs  <- rle(passable)
      r_end <- cumsum(runs$lengths)
      r_str <- r_end - runs$lengths + 1L

      for (k in which(runs$values & runs$lengths >= min_snps)) {
        s         <- r_str[k]
        e         <- r_end[k]
        length_kb <- (bpvec[e] - bpvec[s]) / 1000

        if (length_kb < min_length_kb) next

        seg <- gvec[s:e]
        roh_list[[length(roh_list) + 1L]] <- data.frame(
          iid       = iids[ind_i],
          chr       = ch,
          start_bp  = bpvec[s],
          end_bp    = bpvec[e],
          length_kb = length_kb,
          n_snps    = e - s + 1L,
          n_het     = sum(!is.na(seg) & seg == 1L),
          n_missing = sum(is.na(seg)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(roh_list) == 0) return(.empty_roh_df())
  do.call(rbind, roh_list)
}


# ---- Internal: return an empty but correctly structured ROH data.frame ----
.empty_roh_df <- function() {
  data.frame(
    iid = character(), chr = character(),
    start_bp = integer(), end_bp = integer(),
    length_kb = numeric(), n_snps = integer(),
    n_het = integer(), n_missing = integer(),
    stringsAsFactors = FALSE
  )
}
