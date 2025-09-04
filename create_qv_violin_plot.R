#!/usr/bin/env Rscript

# Create violin plot visualization for QV distributions
# Shows distributions with colored violin plots

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript create_qv_violin_plot.R <output_prefix> [tsv_files...]\n")
  cat("   or: Rscript create_qv_violin_plot.R <output_prefix> --dir <directory>\n")
  cat("\nExample:\n")
  cat("  Rscript create_qv_violin_plot.R violin_report max_qv_results/*.tsv\n")
  cat("  Rscript create_qv_violin_plot.R violin_report --dir max_qv_results\n")
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

cat("Processing", length(tsv_files), "TSV files for violin plot...\n")

# Function to categorize QV values
categorize_qv <- function(qv) {
  case_when(
    qv <= 17 ~ "Very Low (≤17)",
    qv <= 23 ~ "Low (17-23)",
    qv <= 33 ~ "Mid (23-33)",
    TRUE ~ "High (>33)"
  )
}

# Use grayscale for all violins
get_violin_color <- function(mean_qv) {
  return("gray70")  # Light gray for all violins
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
    n = n(),
    .groups = 'drop'
  )

# Use grayscale for all violins
summary_stats <- summary_stats %>%
  mutate(color = sapply(mean_qv, get_violin_color))

# Sort regions lexicographically
summary_stats <- summary_stats %>%
  arrange(region)

# Merge color info back to all_data
all_data <- all_data %>%
  left_join(summary_stats %>% select(region, color), by = "region")

# Order regions lexicographically
all_data$region <- factor(all_data$region, levels = sort(unique(all_data$region)))

# Create violin plot
p <- ggplot(all_data, aes(x = region, y = avg_qv)) +
  geom_violin(aes(fill = region), 
              alpha = 0.7, 
              scale = "width",
              width = 0.8) +
  geom_boxplot(width = 0.15, 
               fill = "white", 
               alpha = 0.8,
               outlier.size = 1,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = setNames(summary_stats$color, summary_stats$region)) +
  
  # Add horizontal lines for QV thresholds
  geom_hline(yintercept = 17, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_hline(yintercept = 23, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_hline(yintercept = 33, linetype = "dashed", color = "gray60", alpha = 0.5) +
  
  # Add threshold labels
  annotate("text", x = 0.5, y = 17, label = "Very Low/Low", hjust = 0, size = 3, color = "gray40") +
  annotate("text", x = 0.5, y = 23, label = "Low/Mid", hjust = 0, size = 3, color = "gray40") +
  annotate("text", x = 0.5, y = 33, label = "Mid/High", hjust = 0, size = 3, color = "gray40") +
  
  # Labels and theme
  labs(
    title = "Quality Value (QV) Distribution Across Genomic Regions",
    subtitle = paste("Violin plots showing distribution for", nrow(summary_stats), "regions,", 
                    nrow(all_data), "total samples"),
    x = "Genomic Region",
    y = "Maximum Attainable QV"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray80", linetype = "solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 65))

# Add mean values as text annotations
for (i in 1:nrow(summary_stats)) {
  p <- p + annotate("text", 
                   x = i, 
                   y = 62,
                   label = sprintf("%.1f", summary_stats$mean_qv[i]),
                   size = 3,
                   fontface = "bold")
}

# Save plots
pdf_file <- paste0(output_prefix, "_violin_plot.pdf")
cat("\nSaving PDF:", pdf_file, "\n")
ggsave(
  pdf_file,
  plot = p,
  width = max(10, length(unique(all_data$region)) * 0.8),
  height = 8,
  units = "in",
  dpi = 300
)

png_file <- paste0(output_prefix, "_violin_plot.png")
cat("Saving PNG:", png_file, "\n")
ggsave(
  png_file,
  plot = p,
  width = max(10, length(unique(all_data$region)) * 0.8),
  height = 8,
  units = "in",
  dpi = 150,
  bg = "white"
)

# Save summary statistics
summary_file <- paste0(output_prefix, "_violin_summary.csv")
cat("Saving summary:", summary_file, "\n")
write_csv(summary_stats %>% select(-color), summary_file)

# Print summary to console
cat("\n========================================\n")
cat("QV Distribution Summary (Violin Plot)\n")
cat("========================================\n\n")
cat(sprintf("%-30s %6s %8s %8s\n", "Region", "N", "Mean QV", "Median"))
cat(sprintf("%-30s %6s %8s %8s\n", 
            "------------------------------", "------", "--------", "--------"))

for (i in 1:nrow(summary_stats)) {
  cat(sprintf("%-30s %6d %8.1f %8.1f\n",
              summary_stats$region[i],
              summary_stats$n[i],
              summary_stats$mean_qv[i],
              summary_stats$median_qv[i]))
}

cat("\nViolin plots in grayscale (black and white)\n")

cat("\nPlots saved successfully!\n")