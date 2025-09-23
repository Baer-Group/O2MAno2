# dnds_genomewide.R - genome-wide dN/dS for O1 and O2
# Inputs: O1.ANN.tsv, O2.ANN.tsv (SnpEff ANN tables), WS292 CDS FASTA
# Outputs (in dnds_results/):
#   - dnds_counts.csv      (N, S, N/S)
#   - dnds_per_site.csv    (N, S, S_sites, N_sites, dN, dS, dN/dS, CIs)
#   - dnds_barplot.png     (per-site dN/dS with error bars)

suppressWarnings({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
})
cat("Working directory:", getwd(), "\n")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)
})

# -------------------- EDIT PATHS IF NEEDED --------------------
o1_tsv <- "O1.ANN.tsv"
o2_tsv <- "O2.ANN.tsv"
cds_fa <- "c_elegans.PRJNA13758.WS292.CDS_transcripts.fa"
outdir <- "dnds_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# -------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a

# ---- Robust SnpEff TSV reader (header/no-header; multiple ANN cols OK) ----
safe_read_ann <- function(path) {
  stopifnot(file.exists(path))
  first <- readLines(path, n = 1, warn = FALSE)
  has_header <- grepl("^CHROM\\tPOS\\tREF\\tALT\\t", first)
  df <- if (has_header) {
    read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE)
  } else {
    read_tsv(path, col_names = c("CHROM","POS","REF","ALT","ANN"),
             col_types = cols(.default = col_character()), progress = FALSE)
  }
  if (ncol(df) > 5) {                                # collapse multi-ANN columns
    df$ANN <- apply(df[, 5:ncol(df), drop = FALSE], 1, function(v) paste(v, collapse = "\t"))
    df <- df[, 1:5]
  }
  if (!"ANN" %in% names(df)) stop("No ANN column found in: ", path)
  df$POS <- suppressWarnings(as.integer(df$POS))
  df
}

# ---- Parse ANN field and tag cohort ----
parse_ann <- function(df, gen_label) {
  df %>%
    filter(!is.na(ANN), ANN != ".") %>%
    mutate(ANN = gsub("\t", ",", ANN, fixed = TRUE)) %>%  # unify separators
    separate_rows(ANN, sep = ",") %>%
    separate(
      ANN,
      into = c("Allele","Annotation","Impact","Gene","Gene_ID",
               "Feature_Type","Feature_ID","Transcript_BioType",
               "Rank","HGVS_c","HGVS_p","cDNA","CDS","AA","Distance","Errors"),
      sep = "\\|", fill = "right", extra = "drop"
    ) %>%
    mutate(Gen = gen_label)
}

# ---- Classify variants once per locus: N if any nonsyn; else S if any syn ----
classify_calls <- function(ann_df) {
  if (nrow(ann_df) == 0)
    return(tibble(Gen=character(), CHROM=character(), POS=integer(),
                  REF=character(), ALT=character(), label=character()))
  
  ann_df %>%
    # keep coding transcripts only
    dplyr::filter(Feature_Type == "transcript",
                  stringr::str_detect(dplyr::coalesce(Transcript_BioType, ""), "protein_coding")) %>%
    # compute flags *after* filtering so lengths match
    dplyr::mutate(
      is_syn  = stringr::str_detect(dplyr::coalesce(Annotation, ""), "synonymous_variant"),
      is_nons = stringr::str_detect(dplyr::coalesce(Annotation, ""),
                                    "missense_variant|stop_gained|stop_lost|start_lost|start_gained")
    ) %>%
    dplyr::group_by(Gen, CHROM, POS, REF, ALT) %>%
    dplyr::summarise(
      label = dplyr::case_when(
        any(is_nons, na.rm = TRUE) ~ "N",
        any(is_syn,  na.rm = TRUE) ~ "S",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(label))
}


# ---- FASTA reader (supports .fa and .fa.gz) ----
read_fasta <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  seqs <- list(); ids <- character(0)
  cur_id <- NULL; buf <- NULL
  while (length(ln <- readLines(con, n = 1, warn = FALSE))) {
    if (!nzchar(ln)) next
    if (substr(ln, 1, 1) == ">") {
      if (!is.null(cur_id)) {
        seqs[[cur_id]] <- toupper(chartr("U","T", gsub("\\s+", "", paste(buf, collapse = ""))))
      }
      hdr <- substring(ln, 2L)
      tok <- strsplit(hdr, "\\s+")[[1]][1]
      tok <- sub("^(transcript:|Transcript:|CDS:)", "", tok)
      cur_id <- tok; buf <- character(0)
    } else {
      buf <- c(buf, ln)
    }
  }
  if (!is.null(cur_id)) {
    seqs[[cur_id]] <- toupper(chartr("U","T", gsub("\\s+", "", paste(buf, collapse = ""))))
  }
  if (!length(seqs)) stop("No FASTA records in: ", path)
  seqs
}

# ---- Nei-Gojobori site counting (precomputed per codon) ----
.codon_aa <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G"
)

# Precompute S/N contribution for each sense codon (stop codons skipped)
.precompute_SN_tables <- local({
  nucs <- c("A","C","G","T")
  S_tab <- numeric(); N_tab <- numeric()
  for (cod in names(.codon_aa)) {
    aa <- .codon_aa[[cod]]
    if (is.null(aa) || aa == "*") next  # ignore stop codons
    syn <- 0; nons <- 0
    for (pos in 1:3) {
      ref <- substr(cod, pos, pos)
      for (alt in nucs[nucs != ref]) {
        mut <- cod; substr(mut, pos, pos) <- alt
        aa2 <- .codon_aa[[mut]]
        if (is.null(aa2)) next
        if (aa2 == aa) syn <- syn + 1 else nons <- nons + 1
      }
    }
    S_tab[[cod]] <- syn/3; N_tab[[cod]] <- nons/3
  }
  list(S_tab = S_tab, N_tab = N_tab)
})

nei_gojobori_sites_all <- function(fasta_seqs) {
  S <- 0; N <- 0
  S_tab <- .precompute_SN_tables$S_tab
  N_tab <- .precompute_SN_tables$N_tab
  for (seq in fasta_seqs) {
    dna <- toupper(gsub("[^ACGT]", "N", seq))
    L <- (nchar(dna) %/% 3L) * 3L
    if (L < 3L) next
    codons <- substring(dna, seq.int(1L, L, 3L), seq.int(3L, L, 3L))
    keep   <- codons %in% names(S_tab)           # drops ambiguous/stop/invalid codons
    if (any(keep)) {
      S <- S + sum(S_tab[codons[keep]], na.rm = TRUE)
      N <- N + sum(N_tab[codons[keep]], na.rm = TRUE)
    }
  }
  c(S_sites = S, N_sites = N)
}

# ---- Poisson CIs for counts  ----
pois_ci <- function(k) {
  if (is.na(k)) return(c(NA_real_, NA_real_))
  if (k == 0) c(0, qchisq(0.975, 2*(k+1))/2) else c(qchisq(0.025, 2*k)/2, qchisq(0.975, 2*(k+1))/2)
}

# ========================= RUN =========================
message("Reading SnpEff tables.")
o1 <- safe_read_ann(o1_tsv) %>% parse_ann("O1")
o2 <- safe_read_ann(o2_tsv) %>% parse_ann("O2")

ann <- bind_rows(o1, o2)
if (nrow(ann) == 0) stop("No ANN rows parsed. Check that O1.ANN.tsv / O2.ANN.tsv contain ANN entries.")

calls <- classify_calls(ann)
if (nrow(calls) == 0) stop("No coding N/S calls found (protein_coding + transcript).")

# ---- counts-based dN/dS ----
counts <- calls %>%
  count(Gen, label) %>%
  pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
  mutate(dNdS_counts = ifelse(S > 0, N / S, NA_real_)) %>%
  arrange(Gen)

write_csv(counts, file.path(outdir, "dnds_counts.csv"))
print(counts)

# ---- per-site dN/dS using CDS FASTA ----
if (!file.exists(cds_fa)) {
  warning("CDS FASTA not found: ", cds_fa, "\nPer-site dN/dS skipped. Counts-only results are in dnds_counts.csv")
  quit(save = "no")
}

message("Computing per-site (Nei-Gojobori) using CDS: ", cds_fa)
cds <- read_fasta(cds_fa)
sites <- nei_gojobori_sites_all(cds)
S_sites <- sites["S_sites"]; N_sites <- sites["N_sites"]


dnds_rates <- counts %>%
  rowwise() %>%
  mutate(
    dN   = N / N_sites,
    dS   = S / S_sites,
    dNdS = dN / dS,
    N_ci = list(pois_ci(N)),
    S_ci = list(pois_ci(S)),
    dN_ci = list(unlist(N_ci) / N_sites),
    dS_ci = list(unlist(S_ci) / S_sites),
    dNdS_ci_low  = (dN_ci[[1]][1]) / (dS_ci[[1]][2]),
    dNdS_ci_high = (dN_ci[[1]][2]) / (dS_ci[[1]][1])
  ) %>%
  ungroup() %>%
  select(Gen, N, S, dN, dS, dNdS, dNdS_ci_low, dNdS_ci_high)

# Katz log-rate-ratio CI for dN/dS
dnds_rates <- dnds_rates %>%
  mutate(
    # rate ratio (you already have dNdS; recompute to be explicit)
    dNdS = (N / N_sites) / (S / S_sites),
    se_log = sqrt(1 / N + 1 / S),          # exposures cancel for Poisson rates
    dNdS_ci_low  = exp(log(dNdS) - 1.96 * se_log),
    dNdS_ci_high = exp(log(dNdS) + 1.96 * se_log)
  )

write_csv(dnds_rates, file.path(outdir, "dnds_per_site.csv"))
print(dnds_rates)

# ---- simple barplot (per-site dN/dS) ----
p <- ggplot(dnds_rates, aes(Gen, dNdS, fill = Gen)) +
  geom_col(width = 0.6, alpha = 0.95) +
  geom_errorbar(aes(ymin = dNdS_ci_low, ymax = dNdS_ci_high), width = 0.2) +
  scale_fill_manual(values = c(O1 = "#0070C0", O2 = "#E46C0A"), guide = "none") +
  labs(y = "dN/dS (per-site; NG)", x = NULL, title = "Genome-wide dN/dS") +
  theme_minimal(base_size = 16)
ggsave(file.path(outdir, "dnds_barplot.png"), p, width = 6, height = 4, dpi = 300)

message("\nDone.\nFiles:\n  - ", file.path(outdir, "dnds_counts.csv"),
        "\n  - ", file.path(outdir, "dnds_per_site.csv"),
        "\n  - ", file.path(outdir, "dnds_barplot.png"))
