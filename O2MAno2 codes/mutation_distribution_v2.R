suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2)
})
# --- optional: set working directory to this script's folder (RStudio only)
try({ setwd(dirname(getActiveDocumentContext()$path)) }, silent = TRUE)
# --- 1) Rebuild points_df (merge E+G into O1/O2, add bins) ---
library(readr); library(dplyr); library(stringr); library(forcats)
getwd()
bin_size <- 250000L  # use your current value

files <- list(
  snp_O1   = c("E.snp.3x_O1.csv",   "G.snp.3x_O1.csv"),
  snp_O2   = c("E.snp.3x_O2.csv",   "G.snp.3x_O2.csv"),
  indel_O1 = c("E.indel.3x_O1.csv", "G.indel.3x_O1.csv"),
  indel_O2 = c("E.indel.3x_O2.csv", "G.indel.3x_O2.csv")
)

normalize_chr <- function(x){
  x <- toupper(gsub("^chr","",as.character(x), ignore.case=TRUE))
  dplyr::recode(x,"1"="I","2"="II","3"="III","4"="IV","5"="V","X"="X", .default=x)
}
read_positions <- function(path){
  df <- readr::read_csv(path, show_col_types = FALSE)
  nm <- names(df)
  chr_col <- nm[stringr::str_detect(nm, regex("^(chrom(osome)?|chr)$", ignore_case=TRUE))][1]
  pos_col <- nm[stringr::str_detect(nm, regex("^(pos(ition)?|bp)$",      ignore_case=TRUE))][1]
  if (is.na(chr_col) || is.na(pos_col)) stop(sprintf("Missing Chromosome/Position in %s", path))
  dplyr::transmute(df,
                   Chromosome = normalize_chr(.data[[chr_col]]),
                   Position   = as.integer(.data[[pos_col]])
  ) |> dplyr::filter(!is.na(Position))
}
prep_one <- function(paths, type_label, gen_label, bin_size){
  dplyr::bind_rows(lapply(paths, read_positions)) |>
    dplyr::distinct(Chromosome, Position, .keep_all = TRUE) |>
    dplyr::mutate(
      Type = type_label,
      Gen  = gen_label,
      Bin  = (Position %/% bin_size) * bin_size
    )
}

points_df <- dplyr::bind_rows(
  prep_one(files$snp_O1,   "SNV",   "O1", bin_size),
  prep_one(files$snp_O2,   "SNV",   "O2", bin_size),
  prep_one(files$indel_O1, "Indel", "O1", bin_size),
  prep_one(files$indel_O2, "Indel", "O2", bin_size)
) |>
  dplyr::mutate(
    Chromosome = forcats::fct_relevel(Chromosome, c("I","II","III","IV","V","X")),
    Type = factor(Type, levels = c("SNV","Indel")),
    Gen  = factor(Gen,  levels = c("O1","O2"))
  )
library(ggplot2); library(scales)

counts_df_bar <- points_df |>
  dplyr::group_by(Chromosome, Type, Gen, Bin) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop")

bar_width_frac <- 0.80   # 0–1 of bin width
x_tick_mb      <- 5

counts_df_bar$Count_adj <- ifelse(counts_df_bar$Gen == "O2", counts_df_bar$Count * 0.5, counts_df_bar$Count)


p_bar <- ggplot(counts_df_bar, aes(x = Bin, y = Count_adj, fill = Gen)) +
  geom_col(
    width = bar_width_frac * bin_size,
    position = position_dodge2(
      width = bin_size, preserve = "single",
      padding = (1 - bar_width_frac) / 2
    )
  ) +
  facet_grid(rows = vars(Type), cols = vars(Chromosome),
             scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(O1 = "#0070C0", O2 = "#E46C0A"),
                    guide = guide_legend(title = NULL, override.aes = list(alpha = 1))) +
  scale_x_continuous(
    breaks = breaks_width(x_tick_mb * 1e6),
    labels = label_number(scale = 1e-6, big.mark = ",")
  ) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Genomic position (in Mb, bin size 250kb) ", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold"),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 20),
    legend.position = "top",
    legend.text = element_text(size = 20)
  )

print(p_bar)
# ggsave("chromosome_panels_bars_O1_vs_O2.png", p_bar, width = 16, height = 9, dpi = 500)
