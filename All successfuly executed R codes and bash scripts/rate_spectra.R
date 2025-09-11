# ---------- Inputs ----------
ref <- 100286401
cov_col <- 9   # 1-based
gen_col <- 10  # capped at 150

# Files (combined E and G within each generation)
f_E_O1 <- "E.3x_O1_all.csv"
f_E_O2 <- "E.3x_O2_all.csv"
f_G_O1 <- "G.3x_O1_all.csv"
f_G_O2 <- "G.3x_O2_all.csv"

# Column index ranges (1-based, inclusive)
snv_cols   <- 2:7          # 6 base-substitution classes
indel_cols <- 12:21        # 10 indel classes (neg_*, pos_*)

# ---------- Libraries ----------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

cap150 <- function(x) pmin(x, 150)

# Per-sample rate table for selected columns (by 1-based index)
rates_for <- function(df, idxs) {
  cov <- df[[cov_col]]
  gens <- cap150(df[[gen_col]])
  denom <- ref * cov * gens
  sub <- df[ , idxs, drop = FALSE]
  sub <- mutate_all(sub, as.numeric)
  as.data.frame(sweep(sub, 1, denom, "/"))
}

# Stack E and G into one matrix of per-sample rates
stack_gen <- function(E_path, G_path, idxs) {
  E <- read_csv(E_path, show_col_types = FALSE)
  G <- read_csv(G_path, show_col_types = FALSE)
  rbind(rates_for(E, idxs), rates_for(G, idxs))
}

# Mean + SEM per class
mean_sem <- function(mat) {
  m <- colMeans(mat, na.rm = TRUE)
  n <- colSums(!is.na(mat))
  s <- apply(mat, 2, sd, na.rm = TRUE) / sqrt(n)
  tibble(category = colnames(mat), mean = as.numeric(m), sem = as.numeric(s))
}

# ---------- SNV spectrum ----------
snv_O1 <- stack_gen(f_E_O1, f_G_O1, snv_cols)
snv_O2 <- stack_gen(f_E_O2, f_G_O2, snv_cols)

snv_O1_df <- mean_sem(snv_O1) %>% mutate(group = "O1")
snv_O2_df <- mean_sem(snv_O2) %>% mutate(group = "O2")
snv_df <- bind_rows(snv_O1_df, snv_O2_df)

write_csv(snv_df, "SNV_rate_spectrum_O1_vs_O2_with_SEM.csv")

# ---------- INDEL spectrum ----------
indel_O1 <- stack_gen(f_E_O1, f_G_O1, indel_cols)
indel_O2 <- stack_gen(f_E_O2, f_G_O2, indel_cols)

indel_O1_df <- mean_sem(indel_O1) %>% mutate(group = "O1")
indel_O2_df <- mean_sem(indel_O2) %>% mutate(group = "O2")
indel_df <- bind_rows(indel_O1_df, indel_O2_df)

write_csv(indel_df, "INDEL_rate_spectrum_O1_vs_O2_with_SEM.csv")

# ---------- Plot helper ----------
plot_grouped <- function(df, title, out_png, rot = 30) {
  ggplot(df, aes(x = category, y = mean, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  position = position_dodge(width = 0.8), width = 0.25) +
    labs(x = NULL, y = "Rate (per bp per generation)", title = title, fill = NULL) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = rot, hjust = 1),
          plot.title = element_text(face = "bold", size = 18))
  ggsave(out_png, width = 10, height = 5, dpi = 300)
}

# ---------- Make figures ----------
plot_grouped(snv_df,   "Six-base substitution rate spectrum (O1 vs O2)", "SNV_rate_spectrum_O1_vs_O2.png")
plot_grouped(indel_df, "Indel (5 ins + 5 del) rate spectrum (O1 vs O2)", "INDEL_rate_spectrum_O1_vs_O2.png")
