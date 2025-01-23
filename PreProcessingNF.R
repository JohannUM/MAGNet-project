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
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

# Load the necessary libraries -------------------------------------------------

library(here)
library(edgeR)
library(ggplot2)
library(limma)

# Load data --------------------------------------------------------------------

setwd(here("data"))

data_samples <- read.csv("MAGNet_PhenoData.csv", row.names = 1)
data_raw_counts <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)

################################################################################
# PREPROCESSING data_samples                                                   #
################################################################################

# Get NF samples ---------------------------------------------------------------

data_samples <- subset(data_samples, etiology == "NF")

# Compute BMI and filter out BMI > 65 ------------------------------------------

data_samples$BMI <- data_samples$weight / ((data_samples$height / 100)^2)
data_samples <- subset(data_samples, BMI <= 65)

# Remove rows with NA in specified columns -------------------------------------


columns_to_check <- c("Diabetes", "age", "race", "gender", "weight", "height", "Library.Pool", "RIN", "TIN.median", "Hypertension", "afib")
data_samples <- data_samples[complete.cases(data_samples[, columns_to_check]), ]

################################################################################
# DATA NORMALIZATION AND TRANSFORMATION                                        #
################################################################################

data_raw_counts <- data_raw_counts[, colnames(data_raw_counts) %in% rownames(data_samples)]

# Statistical testing for confounding factors ----------------------------------

traits <- c("age", "gender", "weight", "height", "Hypertension", "afib")
test_results <- data.frame(Test = character(), Statistic = numeric(), P.Value = numeric(), Significant = logical(), stringsAsFactors = FALSE)

for (trait in traits) {
    if (is.numeric(data_samples[[trait]])) {
        test <- t.test(data_samples[[trait]] ~ data_samples$Diabetes)
        test_results[trait, ] <- c("t-test", test$statistic, test$p.value, test$p.value < 0.05)
    } else {
        test <- chisq.test(table(data_samples[[trait]], data_samples$Diabetes))
        test_results[trait, ] <- c("chi-squared", test$statistic, test$p.value, test$p.value < 0.05)
    }
}

print(test_results)

# Hypertension is significantly different between the groups

# Convert to CPM with TMM and filter -------------------------------------------

dge <- DGEList(counts = data_raw_counts)

dge <- calcNormFactors(dge)

data_tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

data_tmm_cpm_log <- log2(data_tmm_cpm + 1)

keep_genes <- rowSums(data_tmm_cpm_log > 1) >= (ncol(data_tmm_cpm_log) * 0.5)

data_tmm_cpm_log <- data_tmm_cpm_log[keep_genes, ]

data_tmm_cpm_log <- as.data.frame(data_tmm_cpm_log)

# Remove batch effects ---------------------------------------------------------

design <- model.matrix(~ Library.Pool + RIN + TIN.median, data = data_samples)
data_tmm_cpm_log_corrected <- removeBatchEffect(data_tmm_cpm_log, covariates = design[, -1])

# Check for outliers -----------------------------------------------------------

sample_tree <- hclust(dist(t(data_tmm_cpm_log_corrected)), method = "average")
plot(sample_tree, main = "Sample Clustering to Detect Outliers")

# based on clustering P01421 seems to be an outlier, remove them
data_samples <- data_samples[!rownames(data_samples) %in% c("P01421"), ]
data_tmm_cpm_log_corrected <- data_tmm_cpm_log_corrected[, colnames(data_tmm_cpm_log_corrected) %in% rownames(data_samples)]

# Save data --------------------------------------------------------------------

saveRDS(data_tmm_cpm_log_corrected, "NF/data_NF_tmm_cpm_log.RDS")
write.csv(data_samples, file = "NF/data_samples_NF_Diabetes.csv")

################################################################################
