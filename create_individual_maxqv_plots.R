#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(gridExtra)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_individual_maxqv_plots.R <output_dir> <max_qv_dir>")
}

output_dir <- args[1]
max_qv_dir <- args[2]

# Get all max_qv TSV files
tsv_files <- list.files(max_qv_dir, pattern = "*_max_qv\\.tsv$", full.names = TRUE)

if (length(tsv_files) == 0) {
  stop("No TSV files found!")
}

cat("Creating individual plots for", length(tsv_files), "regions...\n\n")

# Function to categorize QV values
categorize_qv <- function(qv) {
  case_when(
    qv <= 17 ~ "Very Low (≤17)",
    qv <= 23 ~ "Low (17-23)",
    qv <= 33 ~ "Mid (23-33)",
    TRUE ~ "High (>33)"
  )
}

# Colors
colors <- c(
  "High (>33)" = "#2ca02c",
  "Mid (23-33)" = "#ffff00",
  "Low (17-23)" = "#ff7f0e",
  "Very Low (≤17)" = "#d62728"
)

# Store all data for combined plot
all_plots <- list()

# Process each file individually
for (i in seq_along(tsv_files)) {
  tsv_file <- tsv_files[i]
  region_name <- basename(tsv_file) %>%
    str_replace("_max_qv\\.tsv$", "")
  
  cat("Processing:", region_name, "\n")
  
  # Read the data
  df <- read_tsv(tsv_file, show_col_types = FALSE)
  df <- df %>%
    mutate(qv_category = categorize_qv(avg_qv))
  
  # Calculate statistics
  mean_qv <- mean(df$avg_qv, na.rm = TRUE)
  median_qv <- median(df$avg_qv, na.rm = TRUE)
  n_samples <- nrow(df)
  
  # Create category percentages
  category_data <- df %>%
    group_by(qv_category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(percentage = 100 * count / n_samples)
  
  # Ensure all categories present
  all_categories <- c("Very Low (≤17)", "Low (17-23)", "Mid (23-33)", "High (>33)")
  category_complete <- data.frame(qv_category = all_categories) %>%
    left_join(category_data, by = "qv_category") %>%
    mutate(
      count = if_else(is.na(count), 0L, count),
      percentage = if_else(is.na(percentage), 0, percentage)
    )
  
  category_complete$qv_category <- factor(
    category_complete$qv_category,
    levels = all_categories
  )
  
  # Create vertical bar plot
  p <- ggplot(category_complete, aes(x = qv_category, y = percentage, fill = qv_category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", percentage)), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    geom_text(aes(label = sprintf("n=%d", count)), 
              vjust = 1.5, size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20),
                      expand = c(0, 0)) +
    labs(
      title = region_name,
      subtitle = sprintf("N=%d | Mean QV=%.1f | Median=%.1f", 
                        n_samples, mean_qv, median_qv),
      x = "",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Store for combined plot
  all_plots[[i]] <- p
  
  # Save individual plot
  output_file <- file.path(output_dir, paste0(region_name, "_maxqv.png"))
  ggsave(output_file, p, width = 6, height = 5, dpi = 150, bg = "white")
  cat("  Saved:", output_file, "\n\n")
}

# Create combined grid plot
if (length(all_plots) > 0) {
  cat("Creating combined grid plot...\n")
  
  # Arrange in 2x3 grid
  ncol_grid <- min(3, length(all_plots))
  nrow_grid <- ceiling(length(all_plots) / ncol_grid)
  
  combined <- do.call(arrangeGrob, c(all_plots, 
                                     list(ncol = ncol_grid, 
                                          nrow = nrow_grid,
                                          top = textGrob("Maximum Attainable QV Distribution by Region",
                                                        gp = gpar(fontsize = 16, fontface = "bold")))))
  
  # Save combined plot
  combined_file <- file.path(output_dir, "all_regions_maxqv_grid.png")
  ggsave(combined_file, combined, 
         width = 6 * ncol_grid, 
         height = 5 * nrow_grid, 
         dpi = 150, bg = "white")
  cat("  Saved combined grid:", combined_file, "\n")
}

cat("\nAll plots created successfully!\n")
