# =============================================================================
# io_plink.R  —  Read PLINK-format SNP data into the standard geno list
#
# STANDARD OUTPUT FORMAT (all io_*.R files return this):
#
#   list(
#     genotypes = integer matrix  [n_ind rows x n_snp cols]
#                 values: 0 = hom-ref, 1 = het, 2 = hom-alt, NA = missing
#     map       = data.frame(chr, snp_id, cm, bp)
#     fam       = data.frame(fid, iid, pid, mid, sex, pheno, population)
#   )
#
# Why 0/1/2 coding?
#   Each value counts the number of copies of the MINOR (less common) allele.
#   0 = both alleles are the major allele (AA)
#   1 = one of each (Aa) — heterozygous
#   2 = both are the minor allele (aa)
#   NA = genotype could not be called (missing data)
# =============================================================================


# -----------------------------------------------------------------------------
# read_plink()
#
# Reads a TEXT PLINK pair: .ped (genotypes) + .map (SNP positions).
#
# .ped layout:  FamID IndID PaternalID MaternalID Sex Phenotype A1 A2 A1 A2 ...
#               6 header cols, then 2 allele columns per SNP
# .map layout:  chromosome  snp_id  genetic_distance_cM  physical_pos_bp
#
# Args:
#   ped_file  - path to .ped file
#   map_file  - path to .map file (defaults to same prefix as ped_file)
# -----------------------------------------------------------------------------
read_plink <- function(ped_file, map_file = NULL) {

  # --- Derive map path if not given ---
  if (is.null(map_file)) {
    map_file <- sub("\\.ped$", ".map", ped_file, ignore.case = TRUE)
  }

  if (!file.exists(ped_file)) stop("PED file not found: ", ped_file)
  if (!file.exists(map_file)) stop("MAP file not found: ", map_file)

  # --- Read the MAP file ---
  map <- utils::read.table(map_file, header = FALSE, stringsAsFactors = FALSE,
                           col.names = c("chr", "snp_id", "cm", "bp"))
  map$chr <- as.character(map$chr)
  n_snp   <- nrow(map)
  message(sprintf("  Read %d SNPs from %s", n_snp, basename(map_file)))

  # --- Read the PED file ---
  ped <- utils::read.table(ped_file, header = FALSE, stringsAsFactors = FALSE,
                           colClasses = "character")
  n_ind <- nrow(ped)
  message(sprintf("  Read %d individuals from %s", n_ind, basename(ped_file)))

  expected_cols <- 6 + 2 * n_snp
  if (ncol(ped) != expected_cols) {
    stop(sprintf("PED has %d columns; expected %d (6 header + 2 x %d SNPs)",
                 ncol(ped), expected_cols, n_snp))
  }

  # --- Extract family / individual information ---
  fam <- data.frame(
    fid        = ped[, 1],
    iid        = ped[, 2],
    pid        = ped[, 3],
    mid        = ped[, 4],
    sex        = suppressWarnings(as.integer(ped[, 5])),
    pheno      = suppressWarnings(as.numeric(ped[, 6])),
    population = NA_character_,
    stringsAsFactors = FALSE
  )

  # --- Recode allele pairs to 0/1/2 ---
  # Extract the allele block (columns 7 onward)
  allele_mat <- as.matrix(ped[, 7:ncol(ped)])
  # PLINK missing code is "0"; also allow "N", ".", "-"
  allele_mat[allele_mat %in% c("0", "N", ".", "-")] <- NA_character_

  genotypes <- matrix(NA_integer_, nrow = n_ind, ncol = n_snp,
                      dimnames = list(fam$iid, map$snp_id))

  for (j in seq_len(n_snp)) {
    a1 <- allele_mat[, 2*j - 1]
    a2 <- allele_mat[, 2*j]

    # Identify the major allele (most common non-missing allele)
    all_alleles <- c(a1, a2)
    all_alleles <- all_alleles[!is.na(all_alleles)]
    if (length(all_alleles) == 0) next   # monomorphic missing — skip

    # Major allele = most frequent
    freq_table <- sort(table(all_alleles), decreasing = TRUE)
    major      <- names(freq_table)[1]

    # Count minor alleles (0/1/2 per individual)
    g <- integer(n_ind)
    for (i in seq_len(n_ind)) {
      if (is.na(a1[i]) || is.na(a2[i])) {
        g[i] <- NA_integer_
      } else {
        # Each allele that is NOT the major allele counts as 1
        g[i] <- as.integer(a1[i] != major) + as.integer(a2[i] != major)
      }
    }
    genotypes[, j] <- g
  }

  list(genotypes = genotypes, map = map, fam = fam)
}


# -----------------------------------------------------------------------------
# read_plink_binary()
#
# Reads binary PLINK files (.bed + .bim + .fam).
# Requires the BEDMatrix package — much faster for large datasets.
#
# The .bim file is the binary equivalent of .map:
#   chr  snp_id  cm  bp  allele1  allele2
#
# Args:
#   prefix - file path WITHOUT extension (e.g., "data/raw/elk_example")
#            OR the full .bed path (extension stripped automatically)
# -----------------------------------------------------------------------------
read_plink_binary <- function(prefix) {

  # Allow the user to pass either the prefix or the .bed filename
  prefix <- sub("\\.bed$", "", prefix, ignore.case = TRUE)

  for (ext in c(".bed", ".bim", ".fam")) {
    path <- paste0(prefix, ext)
    if (!file.exists(path)) stop("File not found: ", path)
  }

  if (!requireNamespace("BEDMatrix", quietly = TRUE)) {
    stop("The BEDMatrix package is required for binary PLINK files.\n",
         "Install it with:  install.packages('BEDMatrix')")
  }

  bm <- BEDMatrix::BEDMatrix(paste0(prefix, ".bed"))

  bim <- utils::read.table(paste0(prefix, ".bim"), header = FALSE,
                           col.names = c("chr", "snp_id", "cm", "bp", "a1", "a2"),
                           stringsAsFactors = FALSE)
  bim$chr <- as.character(bim$chr)

  fam_raw <- utils::read.table(paste0(prefix, ".fam"), header = FALSE,
                               col.names = c("fid", "iid", "pid", "mid", "sex", "pheno"),
                               stringsAsFactors = FALSE)
  fam_raw$population <- NA_character_

  genotypes <- as.matrix(bm[, , drop = FALSE])
  storage.mode(genotypes) <- "integer"
  rownames(genotypes) <- fam_raw$iid
  colnames(genotypes) <- bim$snp_id

  message(sprintf("  Read %d individuals x %d SNPs from binary PLINK",
                  nrow(genotypes), ncol(genotypes)))

  list(
    genotypes = genotypes,
    map       = bim[, c("chr", "snp_id", "cm", "bp")],
    fam       = fam_raw
  )
}
