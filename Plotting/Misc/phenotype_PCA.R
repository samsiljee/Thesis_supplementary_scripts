# Small script to run PCA anlysis on phenotype data summary

# Packages
library(dplyr) # Data manipulation
library(ggplot2) # Plotting

# Files and directories
input_dir <- paste0(getwd(), "/Misc/data/")
output_dir <- paste0(getwd(), "/produced_plots/PCA/")
phenotype_summary_file <- "phenotype_summary.csv"

# Load data
# Load phenotype summary
phenotype_summary <- read.csv(paste0(input_dir, phenotype_summary_file)) %>%
  dplyr::select(-c(Donor, Condition, Well, Index))

# Log transform the data
phenotype_summary[, names(phenotype_summary) != "Sample"] <- log10(phenotype_summary[, names(phenotype_summary) != "Sample"])

# Get separate annotations from summary
annotations <- read.csv(paste0(input_dir, phenotype_summary_file)) %>%
  dplyr::select(c(Donor, Condition, Well, Sample, Index))

# Load functions and custom colours
source(paste0(input_dir, "../../Functions.R"))

# Transform to numeric matrix for PCA analysis
phenotype_mat <- phenotype_summary %>%
  dplyr::select(-"Sample") %>%
  as.matrix() %>%
  t()

# Add back sample names as column names
colnames(phenotype_mat) <- annotations$Index

# Impute missing values and infinite values with row means
k <- which(is.na(phenotype_mat) | phenotype_mat == Inf | phenotype_mat == -Inf, arr.ind=TRUE)

# Calculate row means excluding both NA and infinite values
row_means <- apply(phenotype_mat, 1, function(x) mean(x[is.finite(x)], na.rm=TRUE))

# Replace problematic values with row means
phenotype_mat[k] <- row_means[k[,1]]

# Run PCA analysis
phenotype_PCA <- prcomp(t(phenotype_mat), center = TRUE, scale. = TRUE)

# PCA plot PC1 vs PC2
merge(
  phenotype_PCA$x,
  annotations,
  by.x = "row.names",
  by.y = "Index"
) %>%
  ggplot(aes(x = PC1, y = PC2, colour = Condition, shape = Donor)) +
  geom_point() +
  stat_ellipse(aes(group = Condition)) +
  labs(
    x = paste0(
      "PC1 (",
      round(summary(phenotype_PCA)$importance[2, 1] * 100, digits = 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(phenotype_PCA)$importance[2, 2] * 100, digits = 1),
      "%)"
    ),
    title = "PCA plot - phenotype metrics",
    color = "Condition",
    shape = "Donor"
  ) +
  theme_bw()

# Save the plot
ggsave(paste0(output_dir, "Phenotype PCA plot.png"), width = 8, height = 5.5)
