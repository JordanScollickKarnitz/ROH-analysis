# =============================================================================
# io_generic.R  —  Read a tab-delimited genotype table into the standard geno list
#
# Use this when you have your genotypes in a spreadsheet-style text file
# already encoded as 0/1/2.
#
# REQUIRED COLUMN ORDER:
#   col 1 : snp_id        (e.g., "ELK_SNP_000001")
#   col 2 : chr           (chromosome; e.g., "1", "X")
#   col 3 : bp            (base-pair position; integer)
#   col 4+ : one column per individual (values must be 0, 1, 2, or NA)
#
# OPTIONAL companion sample file:
#   Must have a column named "iid" matching the genotype column headers.
#   Any other columns (population, sex, etc.) are attached to fam.
# =============================================================================

read_generic <- function(geno_file, sample_file = NULL, sep = "\t") {

  if (!file.exists(geno_file)) stop("Genotype file not found: ", geno_file)

  # Use data.table::fread for speed if available, otherwise read.table
  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- data.table::fread(geno_file, sep = sep, data.table = FALSE,
                            na.strings = c("NA", "", "."))
  } else {
    df <- utils::read.table(geno_file, sep = sep, header = TRUE,
                            stringsAsFactors = FALSE,
                            na.strings = c("NA", "", "."))
  }

  if (ncol(df) < 4) {
    stop("Genotype file must have at least 4 columns: snp_id, chr, bp, + 1 individual")
  }

  # --- Build map from first 3 columns ---
  map <- data.frame(
    chr    = as.character(df[[2]]),
    snp_id = as.character(df[[1]]),
    cm     = 0,
    bp     = as.integer(df[[3]]),
    stringsAsFactors = FALSE
  )

  # --- Extract genotype matrix ---
  ind_cols <- 4:ncol(df)
  genotypes <- as.matrix(df[, ind_cols, drop = FALSE])
  storage.mode(genotypes) <- "integer"
  genotypes <- t(genotypes)   # rows = individuals, cols = SNPs
  rownames(genotypes) <- colnames(df)[ind_cols]
  colnames(genotypes) <- map$snp_id

  # --- Validate genotype values ---
  observed_values <- unique(as.vector(genotypes))
  bad_values <- observed_values[!is.na(observed_values) & !observed_values %in% c(0L, 1L, 2L)]
  if (length(bad_values) > 0) {
    warning("Genotype values outside {0, 1, 2, NA} found: ",
            paste(bad_values, collapse = ", "),
            "\nThese will be treated as missing (NA).")
    genotypes[genotypes %in% bad_values] <- NA_integer_
  }

  # --- Build fam ---
  fam <- data.frame(
    fid        = rownames(genotypes),
    iid        = rownames(genotypes),
    pid        = 0,
    mid        = 0,
    sex        = 0,
    pheno      = -9,
    population = NA_character_,
    stringsAsFactors = FALSE
  )

  # --- Attach sample metadata if supplied ---
  if (!is.null(sample_file)) {
    if (!file.exists(sample_file)) {
      warning("Sample file not found, skipping: ", sample_file)
    } else {
      samp <- utils::read.table(sample_file, header = TRUE, sep = sep,
                                stringsAsFactors = FALSE)
      if (!"iid" %in% colnames(samp)) {
        warning("Sample file has no 'iid' column — metadata not attached.")
      } else {
        fam <- merge(fam[, setdiff(colnames(fam), colnames(samp)[colnames(samp) != "iid"])],
                     samp, by = "iid", all.x = TRUE, sort = FALSE)
      }
    }
  }

  message(sprintf("  Loaded %d individuals x %d SNPs from %s",
                  nrow(genotypes), ncol(genotypes), basename(geno_file)))

  list(genotypes = genotypes, map = map, fam = fam)
}
