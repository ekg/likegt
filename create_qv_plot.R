#!/usr/bin/env Rscript

# Create QV distribution plot from likegt max-qv TSV files
# Produces both PDF and PNG output with stacked bar chart

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript create_qv_plot.R <output_prefix> [tsv_files...]\n")
  cat("   or: Rscript create_qv_plot.R <output_prefix> --dir <directory>\n")
  cat("\nExample:\n")
  cat("  Rscript create_qv_plot.R qv_report max_qv_results/*.tsv\n")
  cat("  Rscript create_qv_plot.R qv_report --dir max_qv_results\n")
  quit(status = 1)
}

output_prefix <- args[1]

# Get list of TSV files
if (length(args) > 2 && args[2] == "--dir") {
  # Read all TSV files from directory
  dir_path <- args[3]
  tsv_files <- list.files(dir_path, pattern = "*_max_qv\\.tsv$", full.names = TRUE)
} else {
  # Use provided file list
  tsv_files <- args[-1]
}

if (length(tsv_files) == 0) {
  stop("No TSV files found!")
}

cat("Processing", length(tsv_files), "TSV files...\n")

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
    str_replace("_max_qv\\.tsv$", "") %>%
    str_replace("\\.tsv$", "")
  
  cat("  Processing:", region_name, "\n")
  
  # Read TSV file
  df <- read_tsv(tsv_file, show_col_types = FALSE)
  
  # Add region and categorize
  df <- df %>%
    mutate(
      region = region_name,
      qv_category = categorize_qv(avg_qv)
    )
  
  all_data <- bind_rows(all_data, df)
}

# Calculate summary statistics by region
summary_stats <- all_data %>%
  group_by(region) %>%
  summarise(
    mean_qv = mean(avg_qv, na.rm = TRUE),
    median_qv = median(avg_qv, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_qv))

# Calculate category percentages for each region
category_data <- all_data %>%
  group_by(region, qv_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(region) %>%
  mutate(
    total = sum(count),
    percentage = 100 * count / total
  ) %>%
  ungroup()

# Ensure all categories are present for all regions
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

# Order regions by mean QV (best to worst)
region_order <- summary_stats$region

# Convert to factor with specific order
category_data_complete$region <- factor(category_data_complete$region, levels = rev(region_order))
category_data_complete$qv_category <- factor(
  category_data_complete$qv_category,
  levels = all_categories
)

# Define colors (green to red gradient)
colors <- c(
  "High (>33)" = "#2ca02c",      # Green
  "Mid (23-33)" = "#ffff00",     # Yellow
  "Low (17-23)" = "#ff7f0e",     # Orange
  "Very Low (≤17)" = "#d62728"   # Red
)

# Create the stacked bar plot (vertical bars)
p <- ggplot(category_data_complete, aes(x = region, y = percentage, fill = qv_category)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = colors,
    name = "Quality Category",
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 105),
    expand = c(0, 0)
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  labs(
    title = "Quality Value (QV) Distribution Across Genomic Regions",
    subtitle = paste("Based on maximum attainable QV from", length(tsv_files), "regions"),
    x = "Genomic Region",
    y = "Percentage of Individuals (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Add mean QV annotations (above bars)
mean_qv_labels <- summary_stats %>%
  mutate(region = factor(region, levels = rev(region_order)))

p <- p + 
  geom_text(
    data = mean_qv_labels,
    aes(x = region, y = 102, label = sprintf("%.1f", mean_qv)),
    inherit.aes = FALSE,
    size = 2.5,
    vjust = 0,
    color = "black",
    fontface = "italic"
  )

# Save as PDF
pdf_file <- paste0(output_prefix, "_qv_distribution.pdf")
cat("\nSaving PDF:", pdf_file, "\n")
ggsave(
  pdf_file,
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

# Save as PNG with white background
png_file <- paste0(output_prefix, "_qv_distribution.png")
cat("Saving PNG:", png_file, "\n")
ggsave(
  png_file,
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  dpi = 150,
  bg = "white"
)

# Create summary statistics table
summary_table <- all_data %>%
  group_by(region) %>%
  summarise(
    n = n(),
    mean_qv = mean(avg_qv, na.rm = TRUE),
    sd_qv = sd(avg_qv, na.rm = TRUE),
    median_qv = median(avg_qv, na.rm = TRUE),
    min_qv = min(avg_qv, na.rm = TRUE),
    max_qv = max(avg_qv, na.rm = TRUE),
    very_low = sum(qv_category == "Very Low (≤17)"),
    low = sum(qv_category == "Low (17-23)"),
    mid = sum(qv_category == "Mid (23-33)"),
    high = sum(qv_category == "High (>33)"),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_qv))

# Save summary statistics
summary_file <- paste0(output_prefix, "_qv_summary.csv")
cat("Saving summary:", summary_file, "\n")
write_csv(summary_table, summary_file)

# Print summary to console
cat("\n========================================\n")
cat("QV Distribution Summary\n")
cat("========================================\n\n")
cat(sprintf("%-30s %6s %8s %8s %8s\n", "Region", "N", "Mean QV", "SD", "Range"))
cat(sprintf("%-30s %6s %8s %8s %8s\n", 
            "------------------------------", "------", "--------", "--------", "--------"))

for (i in 1:nrow(summary_table)) {
  cat(sprintf("%-30s %6d %8.1f %8.1f %4.1f-%-4.1f\n",
              summary_table$region[i],
              summary_table$n[i],
              summary_table$mean_qv[i],
              summary_table$sd_qv[i],
              summary_table$min_qv[i],
              summary_table$max_qv[i]))
}

cat("\n")
cat("Overall statistics:\n")
cat(sprintf("  Total individuals: %d\n", nrow(all_data)))
cat(sprintf("  Overall mean QV: %.1f\n", mean(all_data$avg_qv, na.rm = TRUE)))
cat(sprintf("  Overall SD: %.1f\n", sd(all_data$avg_qv, na.rm = TRUE)))

cat("\nPlots saved successfully!\n")