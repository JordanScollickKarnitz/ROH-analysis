# =============================================================================
# tests/test_roh_detection.R
#
# Basic unit tests — run from the repo root with:
#   Rscript tests/test_roh_detection.R
# =============================================================================

if (basename(getwd()) == "tests") setwd("..")

source("R/roh_detection.R")
source("R/roh_metrics.R")
source("R/pop_compare.R")

pass <- 0
fail <- 0

.test <- function(name, expr) {
  result <- tryCatch(expr, error = function(e) {
    cat(sprintf("  [FAIL] %s\n         Error: %s\n", name, e$message))
    fail <<- fail + 1
    return(FALSE)
  })
  if (isTRUE(result)) {
    cat(sprintf("  [PASS] %s\n", name))
    pass <<- pass + 1
  }
}

cat("Running ROH toolkit tests...\n\n")

# ---- Build toy dataset: 3 individuals, 1 chromosome, 300 SNPs ----
n_snp <- 300
map <- data.frame(
  chr    = "1",
  snp_id = paste0("s", seq_len(n_snp)),
  cm     = 0,
  bp     = seq(1, by = 1e4, length.out = n_snp),
  stringsAsFactors = FALSE
)

# IND1: big ROH from SNP 100 to 250 (all homozygous)
g1 <- rep(1L, n_snp)
g1[100:250] <- 2L

# IND2: no ROH (all het)
g2 <- rep(1L, n_snp)

# IND3: small ROH only (too short to pass min_snps = 80)
g3 <- rep(1L, n_snp)
g3[10:50] <- 0L   # 41 SNPs — below threshold

G <- rbind(g1, g2, g3)
rownames(G) <- c("IND1", "IND2", "IND3")
fam <- data.frame(fid = c("A","A","B"), iid = c("IND1","IND2","IND3"),
                   pid = 0, mid = 0, sex = 0, pheno = -9,
                   population = c("PopA","PopA","PopB"),
                   stringsAsFactors = FALSE)
geno <- list(genotypes = G, map = map, fam = fam)

# ---- roh_detection tests ----
roh_c <- detect_roh_consecutive(geno, min_snps = 80, min_length_kb = 0,
                                  verbose = FALSE)

.test("consecutive: finds exactly 1 ROH",
      nrow(roh_c) == 1)

.test("consecutive: ROH belongs to IND1",
      roh_c$iid == "IND1")

.test("consecutive: ROH has >= 150 SNPs (100 to 250)",
      roh_c$n_snps >= 150)

.test("consecutive: IND2 has no ROH",
      !any(roh_c$iid == "IND2"))

.test("consecutive: short ROH in IND3 filtered out",
      !any(roh_c$iid == "IND3"))

roh_s <- detect_roh_sliding(geno, window_snps = 20, min_snps = 80,
                              min_length_kb = 0, max_het_in_window = 0,
                              max_missing_in_window = 0, verbose = FALSE)

.test("sliding: finds at least 1 ROH for IND1",
      any(roh_s$iid == "IND1"))

.test("sliding: IND2 still has no ROH",
      !any(roh_s$iid == "IND2"))

# ---- roh_metrics tests ----
genome_kb <- (n_snp * 1e4) / 1000   # 300 SNPs * 10 kb = 3000 kb

f <- compute_f_roh(roh_c, genome_length_kb = genome_kb,
                    all_individuals = c("IND1", "IND2", "IND3"))

.test("compute_f_roh: returns 3 rows (all individuals)",
      nrow(f) == 3)

.test("compute_f_roh: IND1 F_ROH > 0",
      f$f_roh[f$iid == "IND1"] > 0)

.test("compute_f_roh: IND2 F_ROH == 0",
      f$f_roh[f$iid == "IND2"] == 0)

.test("compute_f_roh: IND1 ranked highest",
      f$iid[1] == "IND1")

# Inject a larger ROH so it falls into the "medium" length class
roh_c$length_kb <- 2500
by_cl <- summarize_roh_by_length(roh_c)

.test("summarize_roh_by_length: returns expected columns",
      all(c("iid", "length_class", "n", "total_kb") %in% colnames(by_cl)))

.test("summarize_roh_by_length: IND1 in medium class",
      any(by_cl$length_class == "medium" & by_cl$iid == "IND1"))

# ---- pop_compare tests ----
f2 <- compute_f_roh(roh_c, genome_length_kb = genome_kb,
                     all_individuals = c("IND1","IND2","IND3"))
cmp <- compare_f_roh(f2, fam)

.test("compare_f_roh: per_pop has expected columns",
      all(c("population","n","mean_f_roh") %in% colnames(cmp$per_pop)))

.test("compare_f_roh: returns a test result",
      !is.null(cmp$test))

# ---- Summary ----
cat(sprintf("\n%d tests passed, %d failed\n", pass, fail))
if (fail > 0) {
  cat("  Some tests FAILED — check output above.\n")
  quit(status = 1)
} else {
  cat("  All tests PASSED.\n")
}
