suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(forcats)
})

# ---- INPUTS (edit as needed)
csv_o1 <- c("E.snp.3x_O1.csv", "G.snp.3x_O1.csv")
csv_o2 <- c("E.snp.3x_O2.csv", "G.snp.3x_O2.csv")

# If your CSVs have headers exactly: Chromosome, Position, REF, ALT you’re good.
# Otherwise, set the col_names vector accordingly.

read_snp_csv <- function(path, col_names = c("Chromosome","Position","REF","ALT")) {
  df <- read_csv(path, show_col_types = FALSE, col_names = TRUE)
  stopifnot(all(col_names %in% names(df)))
  df %>%
    transmute(
      Chromosome = as.character(.data[[col_names[1]]]),
      Position   = as.integer(.data[[col_names[2]]]),
      REF        = as.character(.data[[col_names[3]]]),
      ALT        = as.character(.data[[col_names[4]]])
    )
}

normalize_chr <- function(x) {
  x <- toupper(gsub("^CHR", "", x))
  recode(x, "1"="I","2"="II","3"="III","4"="IV","5"="V", .default = x)
}

prep_gen <- function(paths, gen) {
  bind_rows(lapply(paths, read_snp_csv)) %>%
    mutate(
      Chromosome = normalize_chr(Chromosome),
      REF = toupper(REF), ALT = toupper(ALT)
    ) %>%
    # keep only SNVs (single base ref/alt; drop indels/multiallelic here)
    filter(str_detect(REF, "^[ACGT]$"),
           str_detect(ALT, "^[ACGT]$")) %>%
    distinct(Chromosome, Position, REF, ALT, .keep_all = TRUE) %>%
    mutate(Gen = gen)
}

o1 <- prep_gen(csv_o1, "O1")
o2 <- prep_gen(csv_o2, "O2")

# write minimal VCFs for SnpEff (no samples; just sites)
write_min_vcf <- function(df, outfile) {
  df <- df %>% arrange(Chromosome, Position, REF, ALT)
  hdr <- c(
    "##fileformat=VCFv4.2",
    '##INFO=<ID=SRC,Number=1,Type=String,Description="source (O1/O2)">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )
  con <- gzfile(outfile, "w")
  writeLines(hdr, con)
  if (nrow(df) > 0) {
    lines <- sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\tSRC=%s",
                     df$Chromosome, df$Position, df$REF, df$ALT, df$Gen)
    writeLines(lines, con)
  }
  close(con)
  message("Wrote: ", outfile)
}

write_min_vcf(o1, "O1.snps.min.vcf.gz")
write_min_vcf(o2, "O2.snps.min.vcf.gz")
