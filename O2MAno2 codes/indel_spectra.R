# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
library(scales)
library(matrixStats)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
# Files (E + G, O1 + O2)
f <- list(
  E_O1 = "E.indel.3x_O1_total.csv",
  G_O1 = "G.indel.3x_O1_total.csv",
  E_O2 = "E.indel.3x_O2_total.csv",
  G_O2 = "G.indel.3x_O2_total.csv"
)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Column sets (5 deletion bins, 5 insertion bins)
del_bins <- c("neg_1","neg_2","neg_5_3","neg_10_6","neg_inf_10")
ins_bins <- c("pos_1","pos_2","pos_3_5","pos_6_10","pos_10_inf")
bin_labels <- c("1 bp","2 bp","3-5 bp","6-10 bp","10+ bp")


read_tbl <- function(path) readr::read_csv(path, show_col_types = FALSE)

# bind E and G within each generation (row-wise = samples)
O1 <- bind_rows(read_tbl(f$E_O1), read_tbl(f$G_O1))
O2 <- bind_rows(read_tbl(f$E_O2), read_tbl(f$G_O2))

# Turn counts into per-sample proportions within the 5 bins
row_prop <- function(df, bins) {
  num <- df[, bins]
  tot <- rowSums(num, na.rm = TRUE)
  sweep(num, 1, ifelse(tot == 0, NA, tot), "/")
}

O1_del_prop <- row_prop(O1, del_bins)
O2_del_prop <- row_prop(O2, del_bins)
O1_ins_prop <- row_prop(O1, ins_bins)
O2_ins_prop <- row_prop(O2, ins_bins)

# mean + SEM across samples for each bin
mean_sem <- function(prop_df) {
  m <- colMeans(prop_df, na.rm = TRUE)
  n <- colSums(!is.na(prop_df))
  s <- apply(prop_df, 2, sd, na.rm = TRUE) / sqrt(n)
  tibble(bin = bin_labels, mean = as.numeric(m), sem = as.numeric(s))
}

del_O1 <- mean_sem(O1_del_prop) %>% mutate(group = "O1")
del_O2 <- mean_sem(O2_del_prop) %>% mutate(group = "O2")
ins_O1 <- mean_sem(O1_ins_prop) %>% mutate(group = "O1")
ins_O2 <- mean_sem(O2_ins_prop) %>% mutate(group = "O2")

del_df <- bind_rows(del_O1, del_O2)
ins_df <- bind_rows(ins_O1, ins_O2)

# Save tables
write_csv(del_df, "Deletion_proportions_O1_vs_O2_with_SEM.csv")
write_csv(ins_df, "Insertion_proportions_O1_vs_O2_with_SEM.csv")

# Plot helper: grouped bars with error bars
plot_spectrum <- function(df, title, out_png) {
  # --- enforce bin order (watch the dash character!) ---
  order_levels <- c("1 bp","2 bp","3-5 bp","6-10 bp","10+ bp")
  df$bin <- factor(df$bin, levels = order_levels, ordered = TRUE)
  ggplot(df, aes(x = bin, y = mean, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  position = position_dodge(width = 0.8), width = 0.25) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Proportion", title = title, fill = NULL) +
    theme_bw(base_size = 18) +        # <— raises the default font for everything
    theme(
      text        = element_text(size = 18),
      axis.title  = element_text(size = 18),
      axis.text   = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title= element_text(size = 16),
      plot.title  = element_text(size = 22, face = "bold")
    )
  ggsave(out_png, width = 8, height = 4.6, dpi = 400)
}

plot_spectrum(del_df, "Deletion spectrum ( with SEM )",
              "Deletion_Spectrum_with_SEM.png")
plot_spectrum(ins_df, "Insertion spectrum ( with SEM )",
              "Insertion_Spectrum_with_SEM.png")

