#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir <- ifelse(length(args) >= 1, args[[1]], "results/validation_experiments")
out_dir <- ifelse(length(args) >= 2, args[[2]], file.path(in_dir, "figures"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pub_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

# 1) Synteny nulls
synteny_file <- file.path(in_dir, "01_synteny_null.json")
synteny <- fromJSON(synteny_file, simplifyDataFrame = TRUE) %>%
  as_tibble() %>%
  mutate(
    inversion_window = paste0(inv_min, "-", inv_max),
    n_chr = factor(n_chr),
    fusions = factor(fusions)
  )

p1 <- ggplot(
  synteny,
  aes(x = inversions, y = mean_recovery, color = inversion_window, group = inversion_window)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(n_chr ~ fusions, labeller = label_both) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Synteny recovery under simulated inversions/fusions",
    subtitle = "Mean adjacency recovery across parameter combinations",
    x = "Number of inversions",
    y = "Recovered reference adjacencies",
    color = "Inversion size"
  ) +
  pub_theme

ggsave(file.path(out_dir, "figure_synteny_null_recovery.pdf"), p1, width = 10.5, height = 6.8)
ggsave(file.path(out_dir, "figure_synteny_null_recovery.png"), p1, width = 10.5, height = 6.8, dpi = 320)

# 2) Phylogenetic permutation summary
perm_file <- file.path(in_dir, "02_phylo_permutation.json")
perm <- fromJSON(perm_file)
perm_df <- tibble(
  metric = c("Observed\nsocial-solitary effect", "Mean\npermuted effect"),
  value = c(perm$observed_effect, perm$perm_mean_effect)
)

p2 <- ggplot(perm_df, aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.5, size = 3.8) +
  labs(
    title = "Phylogenetic permutation summary",
    subtitle = sprintf("Two-sided empirical P = %.4g (%s permutations)", perm$perm_pvalue_two_sided, format(perm$n_permutations, big.mark=",")),
    x = NULL,
    y = "Effect size"
  ) +
  pub_theme

ggsave(file.path(out_dir, "figure_phylo_permutation_summary.pdf"), p2, width = 7.5, height = 5.5)
ggsave(file.path(out_dir, "figure_phylo_permutation_summary.png"), p2, width = 7.5, height = 5.5, dpi = 320)

# 3) Gene-family simulation contrast
families_file <- file.path(in_dir, "03_gene_family_simulations.json")
families <- fromJSON(families_file)

fam_df <- tibble(
  model = c("Neutral", "Hypothesis"),
  delta_social_minus_solitary = c(
    families$neutral_model$delta_social_minus_solitary,
    families$hypothesis_model$delta_social_minus_solitary
  ),
  loss_modifier = c(
    families$neutral_model$eusocial_loss_modifier,
    families$hypothesis_model$eusocial_loss_modifier
  )
)

p3 <- ggplot(fam_df, aes(x = model, y = delta_social_minus_solitary, fill = model)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("Δ = %.3f", delta_social_minus_solitary)), vjust = -0.45, size = 3.8) +
  labs(
    title = "Gene-family simulations detect modeled sociality effect",
    subtitle = "Reduced eusocial loss rate produces larger social-solitary divergence",
    x = NULL,
    y = "Δ mean family size (social - solitary)"
  ) +
  pub_theme

ggsave(file.path(out_dir, "figure_gene_family_delta.pdf"), p3, width = 7.8, height = 5.5)
ggsave(file.path(out_dir, "figure_gene_family_delta.png"), p3, width = 7.8, height = 5.5, dpi = 320)

# 4) Hamilton-style linkage simulation
ham_file <- file.path(in_dir, "04_hamilton_synteny_simulation.json")
ham <- fromJSON(ham_file, simplifyDataFrame = TRUE) %>% as_tibble()

p4 <- ham %>%
  mutate(margin_sign = if_else(hamilton_margin >= 0, "Positive", "Negative")) %>%
  ggplot(aes(x = recombination_rate, y = final_freq, color = margin_sign)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("margin=%.3f", hamilton_margin)), vjust = -0.8, size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("Positive" = "#1B9E77", "Negative" = "#D95F02")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Hamilton-style simulation across recombination",
    subtitle = "Final social-allele frequency varies with linkage-mediated margin",
    x = "Recombination rate",
    y = "Final social-allele frequency",
    color = "b·R_syn - c"
  ) +
  pub_theme

ggsave(file.path(out_dir, "figure_hamilton_linkage.pdf"), p4, width = 8.6, height = 5.8)
ggsave(file.path(out_dir, "figure_hamilton_linkage.png"), p4, width = 8.6, height = 5.8, dpi = 320)

# Export a compact table for manuscript text/supplement.
summary_tbl <- tibble(
  metric = c(
    "Permutation observed effect",
    "Permutation empirical p",
    "Gene-family Δ neutral",
    "Gene-family Δ hypothesis"
  ),
  value = c(
    perm$observed_effect,
    perm$perm_pvalue_two_sided,
    families$neutral_model$delta_social_minus_solitary,
    families$hypothesis_model$delta_social_minus_solitary
  )
)

write_csv(summary_tbl, file.path(out_dir, "validation_figure_summary_table.csv"))

message("Publication-style figures written to: ", out_dir)
