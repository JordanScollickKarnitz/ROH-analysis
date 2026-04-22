# =============================================================================
# io_vcf.R  —  Read VCF files into the standard geno list
#
# Requires the vcfR package (CRAN):  install.packages("vcfR")
#
# VCF (Variant Call Format) is a text file where each row is a SNP and
# each column after the FORMAT column is an individual's genotype.
# Genotypes are written like "0/0" (homozygous ref), "0/1" (het),
# "1/1" (homozygous alt). Phased genotypes use "|" instead of "/".
#
# This reader:
#   - Drops multi-allelic sites (more than one ALT allele) by default
#   - Converts "./." and "." to NA (missing)
#   - Recodes as dosage of the ALT allele (0/1/2)
# =============================================================================

read_vcf <- function(vcf_file, biallelic_only = TRUE, verbose = TRUE) {

  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("The vcfR package is required for reading VCF files.\n",
         "Install it with:  install.packages('vcfR')")
  }
  if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)

  if (verbose) message("Reading VCF: ", basename(vcf_file))
  vcf <- vcfR::read.vcfR(vcf_file, verbose = verbose)

  # The @fix slot contains the fixed fields: CHROM, POS, ID, REF, ALT, ...
  fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

  # --- Drop multi-allelic SNPs ---
  if (biallelic_only) {
    is_multi <- grepl(",", fix$ALT, fixed = TRUE)
    if (sum(is_multi) > 0) {
      message(sprintf("  Dropping %d multi-allelic sites", sum(is_multi)))
      vcf <- vcf[!is_multi, ]
      fix <- fix[!is_multi, ]
    }
  }

  # --- Extract genotype strings and convert to 0/1/2 ---
  # vcfR::extract.gt returns a matrix of STRINGS like "0/0", "0/1", "1/1"
  # We compute dosage = number of ALT alleles
  gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)

  # Function to convert one genotype string to integer dosage
  gt_to_dosage <- function(x) {
    if (is.na(x) || x %in% c("./.", ".|.", ".")) return(NA_integer_)
    parts <- strsplit(x, "[/|]", perl = TRUE)[[1]]
    if (any(parts == ".")) return(NA_integer_)
    sum(as.integer(parts))
  }

  # Apply over entire matrix (SNPs as rows, samples as cols in vcfR)
  dose <- apply(gt, c(1, 2), gt_to_dosage)
  dose <- t(dose)   # transpose: now individuals are rows, SNPs are cols
  storage.mode(dose) <- "integer"

  # --- Build map ----
  map <- data.frame(
    chr    = as.character(fix$CHROM),
    snp_id = ifelse(is.na(fix$ID) | fix$ID == ".",
                    paste0(fix$CHROM, ":", fix$POS),
                    fix$ID),
    cm     = 0,
    bp     = as.integer(fix$POS),
    stringsAsFactors = FALSE
  )

  # --- Build fam ---
  fam <- data.frame(
    fid        = colnames(gt),
    iid        = colnames(gt),
    pid        = 0,
    mid        = 0,
    sex        = 0,
    pheno      = -9,
    population = NA_character_,
    stringsAsFactors = FALSE
  )

  rownames(dose) <- fam$iid
  colnames(dose) <- map$snp_id

  if (verbose) message(sprintf("  %d individuals x %d SNPs retained",
                               nrow(dose), ncol(dose)))

  list(genotypes = dose, map = map, fam = fam)
}
