# =============================================================================
# simulate_elk_data.R
#
# Generates a synthetic SNP dataset based on North American elk (Cervus canadensis).
#
# Three populations with distinct inbreeding signatures:
#
#   ROCKY  - Rocky Mountain elk (Cervus canadensis nelsoni)
#            Large, panmictic population — baseline diversity, few ROH
#
#   TULE   - Tule elk (Cervus canadensis nannodes)
#            California endemic; extreme historical bottleneck (~1880s near
#            extinction to ~22 individuals). Many MEDIUM ROH from ancient
#            bottleneck. Recovering but still small Ne.
#
#   RANCH  - Ranch/captive elk
#            Small enclosed herds with managed (sometimes accidental)
#            close-kin mating. Many LONG ROH from recent inbreeding.
#
# Output files (in data/raw/):
#   elk_example.ped         PLINK text genotypes
#   elk_example.map         SNP map
#   elk_samples.tsv         Sample metadata
# =============================================================================

set.seed(2024)

# ---- Parameters ----
n_chr        <- 6          # elk have 2n = 68; we simulate 6 autosomes
snps_per_chr <- 2500       # 2500 SNPs / chr
snp_spacing  <- 4e4        # 1 SNP per 40 kb -> ~100 Mb chromosomes
n_snp        <- n_chr * snps_per_chr

# ---- Population sizes ----
n_rocky <- 35
n_tule  <- 30
n_ranch <- 20
pops    <- c(rep("ROCKY", n_rocky), rep("TULE", n_tule), rep("RANCH", n_ranch))
n_ind   <- length(pops)

# ---- Minor allele frequencies ----
# Rocky = most diverse; Tule reduced; Ranch similar to Rocky (recent founder)
maf_rocky <- rbeta(n_snp, 0.6, 2.0)        # broad MAF distribution
maf_tule  <- pmin(maf_rocky * 0.65, 0.50)  # reduced from bottleneck
maf_ranch <- maf_rocky * 0.90               # slight reduction

# ---- Draw genotypes (Hardy-Weinberg) ----
G <- matrix(NA_integer_, nrow = n_ind, ncol = n_snp)
for (i in seq_len(n_ind)) {
  maf <- switch(pops[i],
    ROCKY = maf_rocky,
    TULE  = maf_tule,
    RANCH = maf_ranch
  )
  a1 <- rbinom(n_snp, 1, maf)
  a2 <- rbinom(n_snp, 1, maf)
  G[i, ] <- a1 + a2
}

# ---- Inject MEDIUM ROH into Tule elk (ancient bottleneck signature) ----
# ~3-8 segments of 2-6 Mb per individual, scattered across chromosomes
tule_idx <- which(pops == "TULE")
for (i in tule_idx) {
  n_seg <- sample(4:8, 1)
  for (s in seq_len(n_seg)) {
    ch  <- sample(seq_len(n_chr), 1)
    seg_len_snps <- sample(50:150, 1)   # ~2-6 Mb at 40 kb/SNP
    offset <- (ch - 1) * snps_per_chr
    start  <- offset + sample(1:(snps_per_chr - seg_len_snps), 1)
    allele <- rbinom(1, 1, maf_tule[start]) * 2L
    G[i, start:(start + seg_len_snps - 1)] <- allele
  }
}

# ---- Inject LONG ROH into Ranch elk (recent inbreeding) ----
# 1-2 very long segments (>10 Mb) + 2-4 medium ones per individual
ranch_idx <- which(pops == "RANCH")
for (i in ranch_idx) {
  # Long ROH: 250-400 SNPs = ~10-16 Mb
  n_long <- sample(1:2, 1)
  for (s in seq_len(n_long)) {
    ch  <- sample(seq_len(n_chr), 1)
    seg_len_snps <- sample(250:400, 1)
    offset <- (ch - 1) * snps_per_chr
    start  <- offset + sample(1:(snps_per_chr - seg_len_snps), 1)
    allele <- rbinom(1, 1, 0.3) * 2L
    G[i, start:(start + seg_len_snps - 1)] <- allele
  }
  # Medium ROH on top
  for (s in 1:3) {
    ch  <- sample(seq_len(n_chr), 1)
    seg_len_snps <- sample(60:120, 1)
    offset <- (ch - 1) * snps_per_chr
    start  <- offset + sample(1:(snps_per_chr - seg_len_snps), 1)
    allele <- rbinom(1, 1, 0.3) * 2L
    G[i, start:(start + seg_len_snps - 1)] <- allele
  }
}

# ---- Add a shared ROH island on chr 2 (putative sweep locus) ----
# Present at high frequency in Tule + Ranch (all went through bottleneck)
island_snps <- ((1 * snps_per_chr) + 400):((1 * snps_per_chr) + 550)
island_inds <- c(tule_idx, ranch_idx)
for (i in island_inds) {
  G[i, island_snps] <- 0L  # fix homozygous ref
}

# ---- Light random missingness (~1%) ----
miss_mask <- matrix(rbinom(n_ind * n_snp, 1, 0.01), nrow = n_ind)
G[miss_mask == 1] <- NA_integer_

# ---- Write .map ----
map_df <- data.frame(
  chr    = rep(seq_len(n_chr), each = snps_per_chr),
  snp_id = paste0("ELK_SNP_", sprintf("%06d", seq_len(n_snp))),
  cm     = 0,
  bp     = rep(seq(1, by = snp_spacing, length.out = snps_per_chr), n_chr)
)
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
write.table(map_df, "data/raw/elk_example.map",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Wrote elk_example.map")

# ---- Write .ped ----
iids <- sprintf("ELK_%03d", seq_len(n_ind))
ped_alleles <- matrix("A", nrow = n_ind, ncol = 2 * n_snp)

for (j in seq_len(n_snp)) {
  for (i in seq_len(n_ind)) {
    g <- G[i, j]
    if (is.na(g)) {
      ped_alleles[i, 2*j-1] <- "0"; ped_alleles[i, 2*j] <- "0"
    } else if (g == 0L) {
      ped_alleles[i, 2*j-1] <- "A"; ped_alleles[i, 2*j] <- "A"
    } else if (g == 1L) {
      ped_alleles[i, 2*j-1] <- "A"; ped_alleles[i, 2*j] <- "C"
    } else {
      ped_alleles[i, 2*j-1] <- "C"; ped_alleles[i, 2*j] <- "C"
    }
  }
}

ped_df <- cbind(
  data.frame(fid = iids, iid = iids, pid = 0, mid = 0, sex = 0, pheno = -9),
  as.data.frame(ped_alleles, stringsAsFactors = FALSE)
)
write.table(ped_df, "data/raw/elk_example.ped",
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Wrote elk_example.ped")

# ---- Write sample metadata ----
regions <- c(
  ROCKY = "Western US (Rocky Mountains)",
  TULE  = "California Central Valley",
  RANCH = "Private ranchland (captive)"
)
samples_df <- data.frame(
  iid        = iids,
  population = pops,
  subspecies = ifelse(pops == "TULE", "C. c. nannodes", "C. c. nelsoni"),
  region     = regions[pops],
  sex        = sample(c("M", "F"), n_ind, replace = TRUE),
  year_sampled = sample(2018:2023, n_ind, replace = TRUE),
  stringsAsFactors = FALSE
)
write.table(samples_df, "data/raw/elk_samples.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote elk_samples.tsv")

message(sprintf("\nSimulation complete: %d individuals, %d SNPs, %d chromosomes",
                n_ind, n_snp, n_chr))
message("Expected pattern: RANCH > TULE > ROCKY for F_ROH")
message("  RANCH: many long ROH (recent inbreeding)")
message("  TULE:  many medium ROH (ancient bottleneck)")
message("  ROCKY: few short ROH (large panmictic population)")
