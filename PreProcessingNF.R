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
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")

# Load the necessary libraries -------------------------------------------------

library(here)
library(MatchIt)
library(edgeR)
library(ggplot2)

# Load data --------------------------------------------------------------------

setwd(here("data"))
data_samples <- read.csv("MAGNet_PhenoData.csv", row.names = 1)
data_raw_counts <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)
data_expression_magnet <- readRDS("CPMS_SVA_corrected.RDS")
gene_info_all <- readRDS("gene_info_all.RDS")
gene_lengths <- readRDS("gene_lengths.RDS")

################################################################################
# PREPROCESSING data_samples                                                      #
################################################################################

# Get NF samples ---------------------------------------------------------------

data_samples <- subset(data_samples, etiology == "NF")

# Compute BMI and filter out BMI > 65 ------------------------------------------

data_samples$BMI <- data_samples$weight / ((data_samples$height / 100)^2)
data_samples <- subset(data_samples, BMI <= 65)

################################################################################
# SAMPLE MATCHING                                                              #
################################################################################

# Remove rows with NA in specified columns -------------------------------------

columns_to_check <- c("Diabetes", "age", "race", "gender", "weight", "height", "Library.Pool", "RIN")
data_samples <- data_samples[complete.cases(data_samples[, columns_to_check]), ]

# Convert variables to numeric for matching ------------------------------------

data_samples$Diabetes[which(data_samples$Diabetes == "Yes")] <- 1
data_samples$Diabetes[which(data_samples$Diabetes == "No")] <- 0
data_samples$Diabetes <- as.numeric(as.character(data_samples$Diabetes))

# Optimal matching on a probit PS ----------------------------------------------

m.out <- matchit(
    Diabetes ~ age + race + gender + weight + height,
    data = data_samples,
    method = "optimal",
    distance = "glm",
    link = "probit"
)

# Checking balance after full matching -----------------------------------------

summary(m.out, un = TRUE)
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out,
    type = "density",
    interactive = FALSE,
    which.xs = ~ age + race + gender + weight + height
)
plot(summary(m.out))

# Save matched phenotype data to .csv file -------------------------------------

m.data <- match_data(m.out)

# Restore Diabetes status
m.data$Diabetes[which(m.data$Diabetes == 1)] <- "Yes"
m.data$Diabetes[which(m.data$Diabetes == 0)] <- "No"

data_samples_matched <- as.data.frame(m.data)

################################################################################
# DATA NORMALIZATION AND TRANSFORMATION                                        #
################################################################################

# Statistical testing for confounding factors ----------------------------------

traits <- c("age", "race", "gender", "weight", "height", "Hypertension", "Library.Pool")

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

# Get gene lengths -------------------------------------------------------------

data_raw_counts <- data_raw_counts[, colnames(data_raw_counts) %in% rownames(data_samples_matched)]

dge <- DGEList(counts = data_raw_counts)

dge <- calcNormFactors(dge)

data_tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

data_tmm_cpm_log <- log2(data_tmm_cpm + 1)

keep_genes <- rowSums(data_tmm_cpm_log > 1) >= (ncol(data_tmm_cpm_log) * 0.5)

data_tmm_cpm_log <- data_tmm_cpm_log[keep_genes, ]

data_tmm_cpm_log <- as.data.frame(data_tmm_cpm_log)

sample_tree <- hclust(dist(t(data_tmm_cpm_log)), method = "average")
plot(sample_tree, main = "Sample Clustering to Detect Outliers")

# P01421 P01444 seem to be outliers, remove them
data_samples_matched <- data_samples_matched[!rownames(data_samples_matched) %in% c("P01421", "P01444"), ]
data_tmm_cpm_log <- data_tmm_cpm_log[, colnames(data_tmm_cpm_log) %in% rownames(data_samples_matched)]

saveRDS(data_tmm_cpm_log, "data_NF_tmm_cpm_log.RDS")
write.csv(data_samples_matched, file = "data_samples_NF_Matched_Diabetes.csv")

# design <- model.matrix(~ Hypertension + Library.Pool + RIN, data = data_samples_matched)
# data_tmm_cpm_log_adjusted <- removeBatchEffect(as.data.frame(data_tmm_cpm_log), covariates = design[, -1])
