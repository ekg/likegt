#!/usr/bin/env Rscript

# Create pairwise QV distribution plots from TSV alignment files
# Produces comprehensive visualizations of all pairwise alignments

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(gridExtra)
library(grid)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript create_pairwise_qv_plot.R <output_prefix> [tsv_files...]\n")
  cat("   or: Rscript create_pairwise_qv_plot.R <output_prefix> --pattern <pattern>\n")
  cat("\nExample:\n")
  cat("  Rscript create_pairwise_qv_plot.R pairwise_report pairwise_*_alignments.tsv\n")
  cat("  Rscript create_pairwise_qv_plot.R pairwise_report --pattern 'pairwise_*_alignments.tsv'\n")
  quit(status = 1)
}

output_prefix <- args[1]

# Get list of TSV files
if (length(args) > 2 && args[2] == "--pattern") {
  # Use file pattern
  pattern <- args[3]
  tsv_files <- Sys.glob(pattern)
} else {
  # Use provided file list
  tsv_files <- args[-1]
}

if (length(tsv_files) == 0) {
  stop("No TSV files found!")
}

cat("Processing", length(tsv_files), "pairwise alignment TSV files...\n")

# Function to categorize QV values
categorize_qv <- function(qv) {
  case_when(
    qv <= 17 ~ "Very Low (≤17)",
    qv <= 23 ~ "Low (17-23)",
    qv <= 33 ~ "Mid (23-33)",
    TRUE ~ "High (>33)"
  )
}

# Read and process all TSV files
all_data <- data.frame()

for (tsv_file in tsv_files) {
  # Extract region name from filename
  region_name <- basename(tsv_file) %>%
    str_replace("pairwise_", "") %>%
    str_replace("_pairwise_alignments\\.tsv$", "") %>%
    str_replace("_alignments\\.tsv$", "")
  
  cat("  Processing:", region_name, "\n")
  
  # Read TSV file
  df <- read_tsv(tsv_file, show_col_types = FALSE)
  
  # Filter out self-alignments and add region
  df <- df %>%
    filter(is_self == FALSE | is.na(is_self)) %>%
    mutate(
      region = region_name,
      qv_category = categorize_qv(qv)
    )
  
  all_data <- bind_rows(all_data, df)
}

if (nrow(all_data) == 0) {
  stop("No data found in TSV files!")
}

# Calculate summary statistics by region
summary_stats <- all_data %>%
  group_by(region) %>%
  summarise(
    n_alignments = n(),
    mean_qv = mean(qv, na.rm = TRUE),
    median_qv = median(qv, na.rm = TRUE),
    sd_qv = sd(qv, na.rm = TRUE),
    min_qv = min(qv, na.rm = TRUE),
    max_qv = max(qv, na.rm = TRUE),
    q25 = quantile(qv, 0.25, na.rm = TRUE),
    q75 = quantile(qv, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_qv))

# Define colors (same as max-qv plot)
colors <- c(
  "High (>33)" = "#2ca02c",      # Green
  "Mid (23-33)" = "#ffff00",     # Yellow
  "Low (17-23)" = "#ff7f0e",     # Orange
  "Very Low (≤17)" = "#d62728"   # Red
)

# Create comprehensive plot with 4 panels
fig <- list()

# 1. Violin plot by region
p1 <- ggplot(all_data, aes(x = reorder(region, qv, FUN = median), y = qv, fill = region)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
  coord_flip() +
  labs(
    title = "Pairwise QV Distribution by Region",
    x = "Genomic Region",
    y = "Quality Value (QV)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12)
  )

# 2. Density plot overlay
p2 <- ggplot(all_data, aes(x = qv, color = region)) +
  geom_density(linewidth = 1, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "QV Density Distributions",
    x = "Quality Value (QV)",
    y = "Density",
    color = "Region"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(face = "bold", size = 12)
  ) +
  guides(color = guide_legend(nrow = 2))

# 3. Cumulative distribution
cumulative_data <- all_data %>%
  group_by(region) %>%
  arrange(qv) %>%
  mutate(
    cumulative = seq_along(qv) / n()
  ) %>%
  ungroup()

p3 <- ggplot(cumulative_data, aes(x = qv, y = cumulative, color = region)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Cumulative QV Distribution",
    x = "Quality Value (QV)",
    y = "Cumulative Probability",
    color = "Region"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12)
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(17, 23, 33), linetype = "dotted", alpha = 0.3)

# 4. Category breakdown stacked bar
category_data <- all_data %>%
  group_by(region, qv_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(region) %>%
  mutate(
    total = sum(count),
    percentage = 100 * count / total
  ) %>%
  ungroup()

# Ensure all categories are present
all_categories <- c("Very Low (≤17)", "Low (17-23)", "Mid (23-33)", "High (>33)")
category_data_complete <- expand.grid(
  region = unique(category_data$region),
  qv_category = all_categories,
  stringsAsFactors = FALSE
) %>%
  left_join(category_data, by = c("region", "qv_category")) %>%
  mutate(
    percentage = if_else(is.na(percentage), 0, percentage),
    count = if_else(is.na(count), 0L, count)
  )

category_data_complete$qv_category <- factor(
  category_data_complete$qv_category,
  levels = all_categories
)

p4 <- ggplot(category_data_complete, aes(x = reorder(region, -percentage), y = percentage, fill = qv_category)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = colors,
    name = "QV Category"
  ) +
  labs(
    title = "QV Category Distribution",
    x = "Genomic Region",
    y = "Percentage (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

# Combine all plots
combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2,
                             top = textGrob("Pairwise QV Analysis (All-vs-All Alignments)", 
                                          gp = gpar(fontsize = 16, fontface = "bold")))

# Save combined plot
pdf_file <- paste0(output_prefix, "_pairwise_qv_analysis.pdf")
cat("\nSaving PDF:", pdf_file, "\n")
ggsave(pdf_file, combined_plot, width = 14, height = 10, units = "in", dpi = 300)

png_file <- paste0(output_prefix, "_pairwise_qv_analysis.png")
cat("Saving PNG:", png_file, "\n")
ggsave(png_file, combined_plot, width = 14, height = 10, units = "in", dpi = 150, bg = "white")

# Create overall distribution histogram
p_hist <- ggplot(all_data, aes(x = qv, fill = qv_category)) +
  geom_histogram(bins = 60, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = colors, name = "QV Category") +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  labs(
    title = "Pairwise QV Histogram by Region",
    subtitle = paste("Based on", format(nrow(all_data), big.mark = ","), "total pairwise alignments"),
    x = "Quality Value (QV)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save histogram
hist_pdf <- paste0(output_prefix, "_pairwise_histogram.pdf")
cat("Saving histogram PDF:", hist_pdf, "\n")
ggsave(hist_pdf, p_hist, width = 12, height = 10, units = "in", dpi = 300)

hist_png <- paste0(output_prefix, "_pairwise_histogram.png")
cat("Saving histogram PNG:", hist_png, "\n")
ggsave(hist_png, p_hist, width = 12, height = 10, units = "in", dpi = 150, bg = "white")

# Save summary statistics
summary_file <- paste0(output_prefix, "_pairwise_summary.csv")
cat("Saving summary:", summary_file, "\n")
write_csv(summary_stats, summary_file)

# Print summary to console
cat("\n========================================\n")
cat("Pairwise QV Distribution Summary\n")
cat("========================================\n\n")
cat(sprintf("%-30s %8s %8s %8s %8s\n", "Region", "N_Align", "Mean QV", "Median", "SD"))
cat(sprintf("%-30s %8s %8s %8s %8s\n", 
            "------------------------------", "--------", "--------", "--------", "--------"))

for (i in 1:nrow(summary_stats)) {
  cat(sprintf("%-30s %8d %8.1f %8.1f %8.1f\n",
              summary_stats$region[i],
              summary_stats$n_alignments[i],
              summary_stats$mean_qv[i],
              summary_stats$median_qv[i],
              summary_stats$sd_qv[i]))
}

cat("\nOverall statistics:\n")
cat(sprintf("  Total alignments: %s\n", format(nrow(all_data), big.mark = ",")))
cat(sprintf("  Overall mean QV: %.1f\n", mean(all_data$qv, na.rm = TRUE)))
cat(sprintf("  Overall median QV: %.1f\n", median(all_data$qv, na.rm = TRUE)))
cat(sprintf("  Overall SD: %.1f\n", sd(all_data$qv, na.rm = TRUE)))

# Category breakdown
overall_categories <- all_data %>%
  group_by(qv_category) %>%
  summarise(
    count = n(),
    percentage = 100 * n() / nrow(all_data),
    .groups = 'drop'
  )

cat("\nOverall QV categories:\n")
for (i in 1:nrow(overall_categories)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              overall_categories$qv_category[i],
              overall_categories$count[i],
              overall_categories$percentage[i]))
}

cat("\nPlots saved successfully!\n")