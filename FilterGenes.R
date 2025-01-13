# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

# Load the necessary libraries ------------------------------------------------

library(here)
library(biomaRt)

# Load the data ------------------------------------------------

setwd(here("data"))

# Take preprocessed data here
data_expression <- readRDS("CPMS_SVA_corrected.RDS")
data_samples <- read.csv("phenoData.csv", as.is = T, row.names = 1)
data_exon_lengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, row.names = 1)

data_samples <- data_samples[data_samples$etiology == "DCM", ]
data_expression <- data_expression[, rownames(data_samples)]

# log transform the data
data_log <- log2(data_expression + 1)

if (!all(rownames(data_expression) == rownames(data_exon_lengths))) {
    print("Row names of data_expression and data_exon_lengths do not match or are not in the same order.")
}

# Calculate expression above noise ------------------------------------------------

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    filters = "ensembl_gene_id",
    values = rownames(data_log),
    mart = ensembl
)

rownames(gene_info) <- gene_info$ensembl_gene_id

cpm2fpkm <- function(x) {
    exon_lengths_kb <- data_exon_lengths[, 1] / 1E3
    .t <- 2^(x) / exon_lengths_kb
}
data_expression_fpkm <- cpm2fpkm(data_log)

row_means <- rowMeans(data_expression_fpkm)

genes_on_Y <- rownames(gene_info[gene_info$chromosome_name == "Y", ])

data_expression_fpkm_Y <- data_expression_fpkm[rownames(data_expression_fpkm) %in% genes_on_Y, ]

data_expression_fpkm_Y_female <- data_expression_fpkm_Y[, rownames(data_samples[data_samples$gender == "Female", ])]

noise_mean <- mean(rowMeans(data_expression_fpkm_Y_female))

expression_above_noise <- as.data.frame(t(apply(data_expression_fpkm, 1, function(gene_expression) {
    t_test <- t.test(gene_expression, mu = noise_mean, alternative = "greater")
    above_noise <- t_test$p.value < 0.05
    c(p_value = t_test$p.value, above_noise = above_noise)
})))

data_log_filtered <- data_log[expression_above_noise$above_noise == TRUE, ]

saveRDS(expression_above_noise, "expression_above_noise.RDS")
saveRDS(data_log_filtered, "data_expression_filtered.RDS")
saveRDS(data_samples, "data_samples_filtered.RDS")
