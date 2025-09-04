#!/bin/bash

# Generate all QV plots systematically
# This script creates vertical bar plots for max-qv results
# and optionally runs pairwise analysis

set -e

echo "========================================="
echo "Systematic QV Plot Generation"
echo "========================================="
echo

# Configuration
MAX_QV_DIR="max_qv_results"
PAIRWISE_DIR="pairwise_results"
PLOTS_DIR="qv_plots"
FASTA_DIR="../cositest/alleles_chiara_analysis"
THREADS=8

# Create directories
mkdir -p "$PAIRWISE_DIR"
mkdir -p "$PLOTS_DIR"

# Function to run pairwise analysis if not exists
run_pairwise_if_needed() {
    local fasta="$1"
    local basename=$(basename "$fasta" .fasta.gz)
    local output_file="${PAIRWISE_DIR}/pairwise_${basename}_alignments.tsv"
    
    if [ ! -f "$output_file" ]; then
        echo "  Running pairwise analysis for $basename..."
        python3 analyze_pairwise_qv_simple.py "$fasta" "${PAIRWISE_DIR}/pairwise_${basename}" "$THREADS"
    else
        echo "  Pairwise analysis already exists for $basename"
    fi
}

# Step 1: Check max-qv results
echo "Step 1: Checking max-qv results..."
echo "Found $(ls -1 ${MAX_QV_DIR}/*.tsv 2>/dev/null | wc -l) max-qv result files"
ls -la ${MAX_QV_DIR}/*.tsv 2>/dev/null | tail -10
echo

# Step 2: Generate max-qv vertical bar plot
echo "Step 2: Generating max-qv visualization..."
if ls ${MAX_QV_DIR}/*.tsv 1> /dev/null 2>&1; then
    Rscript create_qv_plot.R "${PLOTS_DIR}/max_qv_report" --dir "$MAX_QV_DIR"
    echo "✅ Max-QV plot saved to ${PLOTS_DIR}/max_qv_report_qv_distribution.png"
else
    echo "⚠️ No max-qv TSV files found in $MAX_QV_DIR"
fi
echo

# Step 3: Run pairwise analysis if needed (optional - comment out if not needed)
echo "Step 3: Checking/running pairwise analyses..."
echo "Note: This uses tree:5:0:0 sparsification for speed"

for fasta in "$FASTA_DIR"/*.fasta.gz; do
    if [ -f "$fasta" ]; then
        basename=$(basename "$fasta" .fasta.gz)
        echo "Checking $basename..."
        
        # Check if pairwise results already exist
        pairwise_file="${PAIRWISE_DIR}/pairwise_${basename}_alignments.tsv"
        if [ ! -f "$pairwise_file" ]; then
            echo "  Running pairwise analysis for $basename..."
            python3 analyze_pairwise_qv_simple.py "$fasta" "${PAIRWISE_DIR}/pairwise_${basename}" "$THREADS"
        else
            echo "  ✓ Pairwise results already exist"
        fi
    fi
done
echo

# Step 4: Generate pairwise visualization if we have data
echo "Step 4: Generating pairwise visualization..."
if ls ${PAIRWISE_DIR}/pairwise_*_alignments.tsv 1> /dev/null 2>&1; then
    Rscript create_pairwise_qv_plot.R "${PLOTS_DIR}/pairwise_report" ${PAIRWISE_DIR}/pairwise_*_alignments.tsv
    echo "✅ Pairwise plot saved to ${PLOTS_DIR}/pairwise_report_pairwise_qv_analysis.png"
else
    echo "⚠️ No pairwise alignment TSV files found"
fi
echo

# Step 5: Create individual region plots
echo "Step 5: Creating individual region comparison plots..."

# Create an R script for individual region plots
cat > create_individual_plots.R << 'EOF'
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
  stop("Usage: Rscript create_individual_plots.R <output_dir> <max_qv_dir>")
}

output_dir <- args[1]
max_qv_dir <- args[2]

# Get all max_qv TSV files
tsv_files <- list.files(max_qv_dir, pattern = "*_max_qv\\.tsv$", full.names = TRUE)

if (length(tsv_files) == 0) {
  stop("No TSV files found!")
}

cat("Creating individual plots for", length(tsv_files), "regions...\n")

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

# Process each file individually
for (tsv_file in tsv_files) {
  region_name <- basename(tsv_file) %>%
    str_replace("_max_qv\\.tsv$", "")
  
  cat("  Processing:", region_name, "\n")
  
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
    geom_text(aes(label = sprintf("%.1f%%", percentage)), 
              vjust = -0.5, size = 3.5) +
    geom_text(aes(label = sprintf("n=%d", count)), 
              vjust = 1.5, size = 3, color = "white") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    labs(
      title = region_name,
      subtitle = sprintf("N=%d samples | Mean QV=%.1f | Median QV=%.1f", 
                        n_samples, mean_qv, median_qv),
      x = "Quality Value Category",
      y = "Percentage of Samples (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.title = element_text(size = 11, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Save individual plot
  output_file <- file.path(output_dir, paste0(region_name, "_qv_distribution.png"))
  ggsave(output_file, p, width = 8, height = 6, dpi = 150, bg = "white")
  cat("    Saved:", output_file, "\n")
}

cat("\nAll individual plots created successfully!\n")
EOF

Rscript create_individual_plots.R "$PLOTS_DIR" "$MAX_QV_DIR"
echo

# Step 6: Summary
echo "========================================="
echo "Summary of Generated Plots"
echo "========================================="
echo
echo "📊 Main visualizations:"
if [ -f "${PLOTS_DIR}/max_qv_report_qv_distribution.png" ]; then
    echo "  ✅ ${PLOTS_DIR}/max_qv_report_qv_distribution.png - Combined max-QV vertical bars"
fi
if [ -f "${PLOTS_DIR}/pairwise_report_pairwise_qv_analysis.png" ]; then
    echo "  ✅ ${PLOTS_DIR}/pairwise_report_pairwise_qv_analysis.png - Pairwise QV analysis"
fi
echo
echo "📈 Individual region plots:"
ls -1 ${PLOTS_DIR}/*_qv_distribution.png 2>/dev/null | grep -v "max_qv_report" | while read file; do
    echo "  ✅ $file"
done
echo
echo "📁 All plots saved in: $PLOTS_DIR/"
echo
echo "Done!"