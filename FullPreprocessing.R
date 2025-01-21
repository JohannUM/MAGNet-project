################################################################################
#                                                                              #
# PRO4002 - Project Period                                                     #
#                                                                              #
# Preprocessing Script                                                         #
#                                                                              #
# K.J. Boessen (i6075128)                                                      #
# N. Imal (i6276013)                                                           #
# J.E. Lottermoser (i6235171)                                                  #
# A. Panghe (i6246854)                                                         #
#                                                                              #
################################################################################

# Install the necessary packages -----------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("MatchIt", quietly = TRUE)) install.packages("MatchIt")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("EDASeq", quietly = TRUE)) BiocManager::install("EDASeq")

# Load the necessary libraries -------------------------------------------------

library(here)
library(MatchIt)
library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)
library(EDASeq)

# Load data --------------------------------------------------------------------

setwd(here("data"))
data_samples <- read.csv("MAGNet_PhenoData.csv", row.names = 1)
data_raw_counts <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)
data_expression_magnet <- readRDS("CPMS_SVA_corrected.RDS")
gene_info_all <- readRDS("gene_info_all.RDS") # Ensembl Backup Data file
gene_lengths <- readRDS("gene_lengths.RDS")

################################################################################
# PREPROCESSING data_samples                                                      #
################################################################################

# Remove Donor, HCM, and PPCM data ---------------------------------------------

data_samples <- subset(data_samples, etiology == "DCM")

# Compute BMI and filter out BMI > 65 ------------------------------------------

data_samples$BMI <- data_samples$weight / ((data_samples$height / 100)^2)
data_samples <- subset(data_samples, BMI <= 65)

################################################################################
# SAMPLE MATCHING                                                              #
################################################################################

# Convert variables to numeric for matching ------------------------------------

# Diabetes
# data_samples_diabetes <- data_samples
# data_samples_diabetes$Diabetes[which(data_samples_diabetes$Diabetes == "Yes")] <- 1
# data_samples_diabetes$Diabetes[which(data_samples_diabetes$Diabetes == "No")] <- 0
# data_samples_diabetes$Diabetes <- as.numeric(as.character(data_samples_diabetes$Diabetes))

# # Ethnicity
# data_samples_ethnicity <- data_samples
# data_samples_ethnicity$race[which(data_samples_ethnicity$race == "AA")] <- 1
# data_samples_ethnicity$race[which(data_samples_ethnicity$race == "Caucasian")] <- 0
# data_samples_ethnicity$race <- as.numeric(as.character(data_samples_ethnicity$race))

# # Optimal matching on a probit PS ----------------------------------------------

# m.out_diabetes <- matchit(
#   Diabetes ~ age + race + gender + weight + height + Hypertension,
#   data = data_samples_diabetes,
#   method = "optimal",
#   distance = "glm",
#   link = "probit"
# )

# m.out_ethnicity <- matchit(
#   race ~ age + gender + weight + height + Diabetes + Hypertension,
#   data = data_samples_ethnicity,
#   method = "optimal",
#   distance = "glm",
#   link = "probit"
# )

# # Checking balance after full matching -----------------------------------------

# summary(m.out_diabetes, un = TRUE)
# plot(m.out_diabetes, type = "jitter", interactive = FALSE)
# plot(m.out_diabetes,
#   type = "density",
#   interactive = FALSE,
#   which.xs = ~ age + race + gender + weight + height + Hypertension
# )
# plot(summary(m.out_diabetes))

# summary(m.out_ethnicity, un = TRUE)
# plot(m.out_ethnicity, type = "jitter", interactive = FALSE)
# plot(m.out_ethnicity,
#   type = "density",
#   interactive = FALSE,
#   which.xs = ~ age + gender + weight + height + Diabetes + Hypertension
# )
# plot(summary(m.out_ethnicity))

# # Save matched phenotype data to .csv file -------------------------------------

# m.data_diabetes <- match_data(m.out_diabetes)
# m.data_ethnicity <- match_data(m.out_ethnicity)

# # Restore Diabetes status
# m.data_diabetes$Diabetes[which(m.data_diabetes$Diabetes == 1)] <- "Yes"
# m.data_diabetes$Diabetes[which(m.data_diabetes$Diabetes == 0)] <- "No"

# # Restore Ethnicity
# m.data_ethnicity$race[which(m.data_ethnicity$race == 1)] <- "AA"
# m.data_ethnicity$race[which(m.data_ethnicity$race == 0)] <- "Caucasian"

# # Save files
# write.csv(m.data_diabetes, file = "MAGNet_data_samples_Matched_Diabetes.csv")
# write.csv(m.data_ethnicity, file = "MAGNet_data_samples_Matched_Ethnicity.csv")

################################################################################
# DATA NORMALIZATION AND TRANSFORMATION                                        #
################################################################################

# Variable of interest switch --------------------------------------------------

# Change argument to data of interest
# data_samples <- as.data.frame(m.data_ethnicity)

# Statistical testing for confounding factors ----------------------------------

write.csv(data_samples, file = "MAGNet_data_samples_Ethnicity.csv")

# Define the traits to test
traits <- c("age", "gender", "weight", "height", "Diabetes", "Hypertension", "Library.Pool")

# Initialize a list to store test results
test_results <- list()

# Perform statistical tests
for (trait in traits) {
  if (is.numeric(data_samples[[trait]])) {
    # Use t-test for numeric traits
    test_results[[trait]] <- t.test(data_samples[[trait]] ~ data_samples$Diabetes)
  } else {
    # Use chi-squared test for categorical traits
    test_results[[trait]] <- chisq.test(table(data_samples[[trait]], data_samples$Diabetes))
  }
}

# Print the test results
for (trait in traits) {
  cat("\nTrait:", trait, "\n")
  print(test_results[[trait]])
}

# Extract protein coding data from Ensembl -------------------------------------

protein_coding <- gene_info_all[gene_info_all$gene_biotype == "protein_coding", ]
num_protein_coding_all <- nrow(protein_coding)
# print(num_protein_coding_all)

num_protein_coding <- nrow(protein_coding[rownames(protein_coding) %in% rownames(data_raw_counts), ])
# print(num_protein_coding)

# Data normalization -----------------------------------------------------------

data_raw_counts <- data_raw_counts[, colnames(data_raw_counts) %in% rownames(data_samples)]
gene_lengths <- gene_lengths[!is.na(as.data.frame(gene_lengths)$length), ]
data_raw_counts <- data_raw_counts[rownames(data_raw_counts) %in% rownames(gene_lengths), ]

# Function to convert raw counts to TPM
counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# Ensure gene lengths are in the same order as raw counts
gene_lengths <- gene_lengths[rownames(data_raw_counts), ]

gene_lengths <- as.data.frame(gene_lengths)

# Convert raw counts to TPM
data_tpm <- apply(data_raw_counts, 2, counts_to_tpm, lengths = gene_lengths$length)

keep_genes_tpm <- rowSums(data_tpm > 1) >= (ncol(data_tpm) * 0.4)

data_tpm <- data_tpm[keep_genes_tpm, ]

data_tpm_log <- log2(data_tpm + 1)

# data_tpm_log_protein <- data_tpm_log[rownames(data_tpm_log) %in% rownames(protein_coding), ]

# saveRDS(data_tpm_log_protein, "data_tpm.RDS")
# all(colnames(data_raw_counts) == rownames(data_samples))

# plot_expression_density(data_tpm_log_protein)

# data_cpm <- cpm(data_raw_counts, log = FALSE)
# keep_genes <- rowSums(data_cpm > 1) >= (ncol(data_cpm) * 0.5)
# keep_genes <- rownames(data_raw_counts) %in% rownames(data_expression_magnet)
# data_raw_counts <- data_raw_counts[keep_genes, ]

# genes <- rownames(data_raw_counts)
# gene_lengths <- getGeneLengthAndGCContent(id = genes, org = "hsa")
# invalid_lengths <- sum(is.na(as.data.frame(gene_lengths)$length))
# print(paste("Number of invalid gene lengths:", invalid_lengths))
# saveRDS(gene_lengths, "gene_lengths.RDS")
# dds <- DESeqDataSetFromMatrix(
#   countData = data_raw_counts,
#   colData = data_samples,
#   design = ~ gender + Library.Pool
# )

# dds <- estimateSizeFactors(dds)

# data_vst <- vst(dds, blind = TRUE)

# data_vst_normalized <- assay(data_vst)

design <- model.matrix(~ age + Library.Pool + Hypertension, data = data_samples)
data_tpm_log_c <- removeBatchEffect(data_tpm_log, covariates = design[, -1])
saveRDS(data_tpm_log_c, "data_tpm_log.RDS")
# plot_expression_density(data_vst_normalized_corrected)

saveRDS(data_vst_normalized_corrected, "final/data_vst_normalized_corrected.RDS")
# CPM + log2 transformation ----------------------------------------------------
# data_cpm <- cpm(data_raw_counts, log = FALSE)

# keep_genes <- rowSums(data_cpm > 1) >= (ncol(data_cpm) * 0.8)
# num_genes_kept <- sum(keep_genes)
# print(num_genes_kept)

# filtered_cpm_matrix <- data_cpm[keep_genes, ]

# data_cpm_magnet <- data_cpm[rownames(data_cpm) %in% rownames(data_expression_magnet), ]

# data_cpm_log_filtered <- log2(filtered_cpm_matrix + 1)
# data_cpm_log_magnet <- log2(data_cpm_magnet + 1)

plot_expression_density <- function(data_cpm_log) {
  expression_values <- as.vector(data_cpm_log)
  density <- ggplot(data.frame(expression_values), aes(x = expression_values)) +
    geom_density(fill = "skyblue", alpha = 0.6) +
    labs(title = "Expression Value Distribution", x = "Log2(CPM)", y = "Density") +
    theme_minimal() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20))

  plot(density)
}

# plot_expression_density(data_cpm_log_filtered)
# plot_expression_density(data_cpm_log_magnet)

sample_tree <- hclust(dist(t(data_tpm_log_protein)), method = "average")
plot(sample_tree, main = "Sample Clustering to Detect Outliers")
# C01997 P01629 C01902 C02660
# data_samples <- data_samples[!rownames(data_samples) %in% c("P01322", "C01902", "C02660"), ]
# data_expression <- data_expression[, colnames(data_expression) %in% rownames(data_samples)]
# P01322 C01902 C02660

# # there is one outlier in the data namely P01629, maybe remove? But matching?

# # correct for confounding factors, maybe add more? Produces negative values because of shifting but is no problem for WGCNA
# design <- model.matrix(~ Library.Pool + gender + Hypertension, data = data_samples)
# adjusted_log_cpm_filtered <- removeBatchEffect(data_cpm_log_filtered, covariates = design[, -1])
# adjusted_log_cpm_magnet <- removeBatchEffect(data_cpm_log_magnet, covariates = design[, -1])

# plot_expression_density(adjusted_log_cpm_filtered)
# plot_expression_density(adjusted_log_cpm_magnet)

# # Save files -------------------------------------------------------------------

# saveRDS(adjusted_log_cpm_filtered, "final/adjusted_log_cpm_filtered.RDS")
# saveRDS(adjusted_log_cpm_magnet, "final/adjusted_log_cpm_magnet.RDS")

################################################################################
