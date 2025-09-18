suppressWarnings({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
})
cat("Working directory:", getwd(), "\n")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(readxl); library(openxlsx); library(purrr)
})

# -------- EDIT THESE --------
in_xlsx  <- "Supplementary_table_S2_List_of_mutations_with_annotation_v3 _2.xlsx"
sheet    <- "Mutations_with_annotation"
out_xlsx <- "Supplementary_table_S2_List_of_mutations_with_annotation_v3.2.xlsx"

ann_files <- c(
  "E.bp.GVCFs.snp.3x.mutation.ann.ANN_expanded.tsv",
  "G.bp.GVCFs.snp.3x.mutation.ann.ANN_expanded.tsv",
  "E.bp.GVCFs.indel.3x.mutation.ann.ANN_expanded.tsv",
  "G.bp.GVCFs.indel.3x.mutation.ann.ANN_expanded.tsv"
)
# ----------------------------

std_chr <- function(x) gsub("^chr","", toupper(as.character(x)))

impact_order <- c("HIGH","MODERATE","LOW","MODIFIER")
effect_order <- c(
  # HIGH
  "frameshift_variant","stop_gained","stop_lost","start_lost","start_gained",
  "transcript_ablation","splice_acceptor_variant","splice_donor_variant",
  # MODERATE
  "splice_region_variant","missense_variant","protein_altering_variant",
  "inframe_insertion","inframe_deletion",
  # LOW
  "synonymous_variant","stop_retained_variant","coding_sequence_variant",
  # UTR
  "5_prime_utr_variant","3_prime_utr_variant",
  # INTRON
  "intron_variant",
  # NEAR-GENE / INTERGENIC
  "upstream_gene_variant","downstream_gene_variant","intergenic_region"
)

read_ann_tsv <- function(path){
  read_tsv(path,
           col_names = c("CHROM","POS","REF","ALT","ANN"),
           show_col_types = FALSE) %>%
    separate(
      ANN,
      into = c("Allele","Annotation","Impact","Gene_name","Gene_id",
               "Feature_type","Feature_id","Transcript_biotype",
               "Rank","HGVS_c","HGVS_p","cDNA","CDS","AA","Distance","Errors"),
      sep = "\\|", fill = "right", extra = "drop"
    )
}

message("Reading ANN_expanded TSVs.")
ann_all <- bind_rows(lapply(ann_files, read_ann_tsv)) %>%
  # keep only rows where ANN allele matches ALT (handles multiallelic)
  filter(is.na(Allele) | Allele == ALT) %>%
  mutate(
    Chromosome = std_chr(CHROM),
    Position   = as.integer(POS),
    Impact_rank = match(Impact, impact_order),
    Impact_rank = ifelse(is.na(Impact_rank), length(impact_order)+1L, Impact_rank),
    Eff_rank    = match(Annotation, effect_order),
    Eff_rank    = ifelse(is.na(Eff_rank), length(effect_order)+1L, Eff_rank)
  ) %>%
  arrange(Chromosome, Position, REF, ALT, Impact_rank, Eff_rank) %>%
  group_by(Chromosome, Position, REF, ALT) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    Chromosome, Position, REF, ALT,
    Annotation_new      = Annotation,
    Putative_impact_new = Impact,
    Gene_name_new       = Gene_name,
    Gene_id_new         = Gene_id
  )

message("New annotations parsed: ", nrow(ann_all))

# --- Load Excel and make the join robust to small header differences ---
wb_df <- read_excel(in_xlsx, sheet = sheet) %>% as_tibble()

# Normalize any hidden spaces in names
names(wb_df) <- names(wb_df) |>
  str_replace_all("\u00A0"," ") |> str_replace_all("\\s+"," ") |> str_trim()

# Helper to pick a column flexibly
pick_col <- function(df, patterns){
  nm <- names(df)
  hits <- which(Reduce(`|`, lapply(patterns, function(p) grepl(p, nm, ignore.case = TRUE, perl = TRUE))))
  if (length(hits)) nm[hits[1]] else NA_character_
}

# Find columns (accept common variants)
ann_name   <- pick_col(wb_df, c("^annotation$","^ann$"))
imp_name   <- pick_col(wb_df, c("^putative[ _]?impact$","^impact$"))
gname_name <- pick_col(wb_df, c("^gene[ _]?name$","^gene$"))
gid_name   <- pick_col(wb_df, c("^gene[ _]?id$","^wb.*id$"))

# Keys must exist as exactly these
if (!all(c("Chromosome","Position","REF","ALT") %in% names(wb_df))){
  stop("Excel must contain columns: Chromosome, Position, REF, ALT")
}

# Keep old versions (no rename, index by the names we found)
wb_df$Annotation_old      <- if (!is.na(ann_name))  wb_df[[ann_name]]  else NA_character_
wb_df$Putative_impact_old <- if (!is.na(imp_name))  wb_df[[imp_name]]  else NA_character_
wb_df$Gene_name_old       <- if (!is.na(gname_name)) wb_df[[gname_name]] else NA_character_
wb_df$Gene_id_old         <- if (!is.na(gid_name))   wb_df[[gid_name]]   else NA_character_

# Standardize keys for join
wb_df <- wb_df %>%
  mutate(Chromosome_key = std_chr(Chromosome),
         Position_key   = as.integer(Position))

# Join and coalesce
out <- wb_df %>%
  left_join(ann_all,
            by = c("Chromosome_key"="Chromosome",
                   "Position_key"  ="Position",
                   "REF","ALT")) %>%
  mutate(
    Annotation      = dplyr::coalesce(.data[["Annotation_new"]],      .data[["Annotation_old"]]),
    Putative_impact = dplyr::coalesce(.data[["Putative_impact_new"]], .data[["Putative_impact_old"]]),
    Gene_name       = dplyr::coalesce(.data[["Gene_name_new"]],       .data[["Gene_name_old"]]),
    Gene_id         = dplyr::coalesce(.data[["Gene_id_new"]],         .data[["Gene_id_old"]])
  )

# QC
qc <- tibble(
  rows_in_excel                 = nrow(wb_df),
  annotation_replaced_from_vcf  = sum(!is.na(out$Annotation_new)),
  rows_without_new_annotation   = sum(is.na(out$Annotation_new))
)

diffs <- out %>%
  filter(!is.na(Annotation_new) & !is.na(Annotation_old) & Annotation_new != Annotation_old) %>%
  select(Chromosome, Position, REF, ALT,
         Annotation_old, Annotation_new,
         Putative_impact_old, Putative_impact_new)

unmatched <- out %>%
  filter(is.na(Annotation_new)) %>%
  select(Chromosome, Position, REF, ALT, Annotation_old, Putative_impact_old)

# Write workbook
wb <- createWorkbook()
addWorksheet(wb, "Updated")
addWorksheet(wb, "Diffs_new_vs_old")
addWorksheet(wb, "Unmatched")
addWorksheet(wb, "QC")

updated <- out %>% select(-ends_with("_key"), -ends_with("_new"), -ends_with("_old"))
writeDataTable(wb, "Updated", updated)
writeDataTable(wb, "Diffs_new_vs_old", diffs)
writeDataTable(wb, "Unmatched", unmatched)
writeDataTable(wb, "QC", qc)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Wrote: ", out_xlsx)


