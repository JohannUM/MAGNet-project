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

# Load the necessary libraries -------------------------------------------------

library(here)
library(edgeR)
library(ggplot2)

# Load data --------------------------------------------------------------------

setwd(here("data"))

data_samples <- read.csv("MAGNet_PhenoData.csv", row.names = 1)
data_raw_counts <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)

################################################################################
# PREPROCESSING data_samples                                                      #
################################################################################

# Remove Donor, HCM, and PPCM data ---------------------------------------------

data_samples <- subset(data_samples, etiology == "DCM")

# Compute BMI and filter out BMI > 65 ------------------------------------------

data_samples$BMI <- data_samples$weight / ((data_samples$height / 100)^2)
data_samples <- subset(data_samples, BMI <= 65)

# Remove rows with NA in specified columns -------------------------------------

columns_to_check <- c("Diabetes", "age", "race", "gender", "weight", "height", "Library.Pool", "RIN")
data_samples <- data_samples[complete.cases(data_samples[, columns_to_check]), ]


################################################################################
# DATA NORMALIZATION AND TRANSFORMATION                                        #
################################################################################

traits <- c("age", "gender", "weight", "height", "Hypertension", "Library.Pool")
test_results <- list()

for (trait in traits) {
    if (is.numeric(data_samples[[trait]])) {
        test_results[[trait]] <- t.test(data_samples[[trait]] ~ data_samples$Diabetes)
    } else {
        test_results[[trait]] <- chisq.test(table(data_samples[[trait]], data_samples$Diabetes))
    }
}

for (trait in traits) {
    cat("\nTrait:", trait, "\n")
    print(test_results[[trait]])
}

data_raw_counts <- data_raw_counts[, colnames(data_raw_counts) %in% rownames(data_samples)]

dge <- DGEList(counts = data_raw_counts)

dge <- calcNormFactors(dge)

data_tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

data_tmm_cpm_log <- log2(data_tmm_cpm + 1)

keep_genes <- rowSums(data_tmm_cpm_log > 1) >= (ncol(data_tmm_cpm_log) * 0.5)

data_tmm_cpm_log <- data_tmm_cpm_log[keep_genes, ]

data_tmm_cpm_log <- as.data.frame(data_tmm_cpm_log)

design <- model.matrix(~ Library.Pool + RIN + Hypertension, data = data_samples)
log_tmm_cpm_corrected <- removeBatchEffect(data_tmm_cpm_log, covariates = design[, -1])

################################################################################
