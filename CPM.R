# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

# Load the necessary libraries ------------------------------------------------

library(here)
library(edgeR)
library(ggplot2)
library(limma)

# Load the data ------------------------------------------------

setwd(here("data"))

data_raw_counts <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)
data_expression <- readRDS("CPMS_SVA_corrected.RDS")
data_samples <- read.csv("MAGNet_PhenoData_Matched.csv", row.names = 1)
gene_info_all <- readRDS("gene_info_all.RDS")

protein_coding <- gene_info_all[gene_info_all$gene_biotype == "protein_coding", ]
num_protein_coding_all <- nrow(protein_coding)
print(num_protein_coding_all)

num_protein_coding <- nrow(protein_coding[rownames(protein_coding) %in% rownames(data_raw_counts), ])
print(num_protein_coding)

data_raw_counts <- data_raw_counts[, colnames(data_raw_counts) %in% rownames(data_samples)]

data_cpm <- cpm(data_raw_counts, log = FALSE)

keep_genes <- rowSums(data_cpm > 1) >= (ncol(data_cpm) / 2)
num_genes_kept <- sum(keep_genes)
print(num_genes_kept)
filtered_cpm_matrix <- data_cpm[keep_genes, ]

magnet_cpm <- data_cpm[rownames(data_cpm) %in% rownames(data_expression), ]

data_cpm_log_filtered <- log2(filtered_cpm_matrix + 1)
data_cpm_log_magnet <- log2(magnet_cpm + 1)

plot_expression_density <- function(data_cpm_log) {
    expression_values <- as.vector(data_cpm_log)
    density <- ggplot(data.frame(expression_values), aes(x = expression_values)) +
        geom_density(fill = "skyblue", alpha = 0.6) +
        labs(title = "Expression Value Distribution", x = "Log2(CPM)", y = "Density") +
        theme_minimal() +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 20))

    plot(density)
}

plot_expression_density(data_cpm_log_filtered)
plot_expression_density(data_cpm_log_magnet)

sample_tree <- hclust(dist(t(data_cpm_log_filtered)), method = "average")
plot(sample_tree, main = "Sample Clustering to Detect Outliers")
# there is one outlier in the data namely P01629, maybe remove? But matching?

# correct for confounding factors, maybe add more? Produces negative values because of shifting but is no problem for WGCNA
design <- model.matrix(~ Library.Pool + gender + Hypertension, data = data_samples)
adjusted_log_cpm_filtered <- removeBatchEffect(data_cpm_log_filtered, covariates = design[, -1])
adjusted_log_cpm_magnet <- removeBatchEffect(data_cpm_log_magnet, covariates = design[, -1])

plot_expression_density(adjusted_log_cpm_filtered)
plot_expression_density(adjusted_log_cpm_magnet)

saveRDS(adjusted_log_cpm_filtered, "final/data_expression_filtered_self.RDS")
saveRDS(adjusted_log_cpm_magnet, "final/data_expression_filtered_magnet.RDS")
