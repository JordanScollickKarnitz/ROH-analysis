# =============================================================================
# run_elk_roh_analysis.R
#
# COMPLETE ELK ROH ANALYSIS — run this top to bottom.
#
# This is the main script students should work through. It:
#   1. Generates simulated elk data (3 populations)
#   2. Loads the data
#   3. Detects ROH using two methods
#   4. Computes F_ROH and length-class summaries
#   5. Compares populations statistically
#   6. Identifies ROH islands
#   7. Produces five figures
#   8. Saves all results to CSV
#
# Time to run: ~1-2 minutes on a laptop
# =============================================================================

cat("====================================================\n")
cat("  ELK RUNS OF HOMOZYGOSITY (ROH) ANALYSIS\n")
cat("  Cervus canadensis — three population comparison\n")
cat("====================================================\n\n")

# ---- 0. Load all functions -------------------------------------------------
source("R/io_plink.R")
source("R/io_vcf.R")
source("R/io_generic.R")
source("R/roh_detection.R")
source("R/roh_metrics.R")
source("R/pop_compare.R")
source("R/plotting.R")

# ---- 1. Generate example data (skip if you already have real data) ----------
cat("STEP 1: Generating simulated elk SNP data...\n")
source("scripts/simulate_elk_data.R")

# ---- 2. Load genotype data -------------------------------------------------
cat("\nSTEP 2: Loading genotype data...\n")
geno <- read_plink("data/raw/elk_example.ped",
                    map_file = "data/raw/elk_example.map")

# Attach population metadata
samples <- read.table("data/raw/elk_samples.tsv", header = TRUE,
                       sep = "\t", stringsAsFactors = FALSE)
geno$fam <- merge(geno$fam[, c("fid","iid","pid","mid","sex","pheno")],
                   samples, by = "iid", all.x = TRUE, sort = FALSE)

# Quick check
cat(sprintf("  Individuals: %d\n", nrow(geno$fam)))
cat(sprintf("  SNPs:        %d\n", ncol(geno$genotypes)))
cat("  Population breakdown:\n")
print(table(geno$fam$population))

# ---- 3. Detect ROH — sliding window ----------------------------------------
cat("\nSTEP 3: Detecting ROH (sliding window method)...\n")
roh <- detect_roh_sliding(
  geno,
  window_snps           = 50,      # window size in SNPs
  min_snps              = 80,      # minimum ROH length in SNPs
  min_length_kb         = 1000,    # minimum ROH length (1 Mb)
  max_het_in_window     = 1,       # allow 1 het per window (error tolerance)
  max_missing_in_window = 2,       # allow 2 missing per window
  verbose               = TRUE
)

cat(sprintf("\n  Total ROH segments found: %d\n", nrow(roh)))
cat(sprintf("  Individuals with at least 1 ROH: %d\n",
            length(unique(roh$iid))))

# Also run the consecutive method for comparison
cat("\nSTEP 3b: Running consecutive method for comparison...\n")
roh_consec <- detect_roh_consecutive(
  geno,
  min_snps      = 80,
  min_length_kb = 1000,
  allow_missing = FALSE,
  verbose       = FALSE
)
cat(sprintf("  Consecutive method found: %d ROH (vs %d from sliding window)\n",
            nrow(roh_consec), nrow(roh)))

# ---- 4. Compute F_ROH -------------------------------------------------------
cat("\nSTEP 4: Computing genomic inbreeding coefficients...\n")

# Elk (Cervus canadensis) autosomal genome ~ 2.5 Gb
# Our simulation uses 6 chromosomes of ~100 Mb each = 600 Mb = 6e5 kb
GENOME_KB <- 6e5   # CHANGE THIS to 2.5e6 for real elk data!

f_roh <- compute_f_roh(roh,
                        genome_length_kb = GENOME_KB,
                        all_individuals  = geno$fam$iid)

by_class <- summarize_roh_by_length(roh)

summary_tbl <- roh_summary_table(roh,
                                  genome_length_kb = GENOME_KB,
                                  fam = geno$fam[, c("iid", "population")])

cat("\nF_ROH summary (top 15 most inbred individuals):\n")
print(head(summary_tbl[, c("iid","population","n_roh","total_roh_kb","f_roh")], 15))

# ---- 5. Compare populations -------------------------------------------------
cat("\n\nSTEP 5: Statistical comparison across populations...\n")
cmp <- compare_f_roh(f_roh, geno$fam)

cat("\nPer-population F_ROH summary:\n")
print(cmp$per_pop)

cat("\nStatistical test result:\n")
print(cmp$test)

if (!is.null(cmp$posthoc)) {
  cat("\nPost-hoc pairwise comparisons (Bonferroni corrected):\n")
  print(cmp$posthoc)
}

cat("\nLength-class comparisons:\n")
lc_cmp <- compare_length_classes(by_class, geno$fam)
print(lc_cmp)

# ---- 6. ROH islands ---------------------------------------------------------
cat("\nSTEP 6: Identifying ROH islands...\n")
islands <- roh_islands(roh,
                        bin_kb        = 250,
                        top_frac      = 0.95,
                        n_individuals = nrow(geno$genotypes))
cat(sprintf("  Found %d candidate ROH islands\n", nrow(islands)))
if (nrow(islands) > 0) {
  cat("\nTop 10 ROH islands:\n")
  print(head(islands, 10))
}

# ---- 7. Save results to CSV ------------------------------------------------
cat("\nSTEP 7: Saving results...\n")
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(roh,          "results/tables/roh_segments.csv",    row.names = FALSE)
write.csv(f_roh,        "results/tables/f_roh.csv",           row.names = FALSE)
write.csv(by_class,     "results/tables/roh_by_length.csv",   row.names = FALSE)
write.csv(summary_tbl,  "results/tables/roh_summary.csv",     row.names = FALSE)
write.csv(cmp$per_pop,  "results/tables/pop_comparison.csv",  row.names = FALSE)
write.csv(islands,      "results/tables/roh_islands.csv",     row.names = FALSE)
cat("  All tables saved to results/tables/\n")

# ---- 8. Figures -------------------------------------------------------------
cat("\nSTEP 8: Generating figures...\n")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("  ggplot2 not installed. Skipping figures.")
  message("  Install with:  install.packages('ggplot2')")
} else {

  # Figure 1: F_ROH boxplot
  p1 <- plot_f_roh_boxplot(f_roh, geno$fam,
                             pop_order = c("ROCKY", "TULE", "RANCH"))
  ggplot2::ggsave("results/figures/fig1_f_roh_boxplot.png",
                   p1, width = 6, height = 5, dpi = 300)

  # Figure 2: ROH length distribution
  p2 <- plot_roh_distribution(roh, fam = geno$fam, by = "population")
  ggplot2::ggsave("results/figures/fig2_roh_length_distribution.png",
                   p2, width = 8, height = 4, dpi = 300)

  # Figure 3: Length-class stacked bars
  p3 <- plot_length_class_stack(by_class, fam = geno$fam)
  ggplot2::ggsave("results/figures/fig3_length_class_stack.png",
                   p3, width = 12, height = 5, dpi = 300)

  # Figure 4: Karyotype for the 6 most inbred individuals
  top6 <- head(f_roh$iid[order(-f_roh$f_roh)], 6)
  p4   <- plot_roh_karyotype(roh, iids = top6,
                               title = "ROH in 6 most inbred elk individuals")
  ggplot2::ggsave("results/figures/fig4_karyotype_top6.png",
                   p4, width = 12, height = 7, dpi = 300)

  # Figure 5: ROH island density on chromosome 2 (injected island location)
  p5 <- plot_roh_island_density(roh, chromosome = "2",
                                 bin_kb = 250,
                                 n_individuals = nrow(geno$genotypes))
  ggplot2::ggsave("results/figures/fig5_roh_island_chr2.png",
                   p5, width = 9, height = 4, dpi = 300)

  cat("  5 figures saved to results/figures/\n")
}

cat("\n====================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("  Check results/tables/ and results/figures/\n")
cat("====================================================\n")
