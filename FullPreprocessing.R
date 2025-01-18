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

# Load the necessary libraries -------------------------------------------------

library(here)
library(MatchIt)
library(edgeR)
library(limma)

################################################################################
# PREPROCESSING PHENODATA                                                      #
################################################################################

# Load data --------------------------------------------------------------------

setwd(here("data"))
phenoData <- read.csv("MAGNet_PhenoData.csv", row.names = 1)

# Remove Donor, HCM, and PPCM data ---------------------------------------------

phenoData <- subset(phenoData, etiology == "DCM")

# Compute BMI and filter out BMI > 65 ------------------------------------------

phenoData$BMI <- phenoData$weight / ((phenoData$height/100)^2)
phenoData <- subset(phenoData, BMI <= 65)

################################################################################
# SAMPLE MATCHING                                                              #
################################################################################

# Convert variables to numeric for matching ------------------------------------

# Diabetes
phenoData1 <- phenoData
phenoData1$Diabetes[which(phenoData1$Diabetes == "Yes")] <- 1
phenoData1$Diabetes[which(phenoData1$Diabetes == "No")] <- 0
phenoData1$Diabetes <- as.numeric(as.character(phenoData1$Diabetes))

# Ethnicity
phenoData2 <- phenoData
phenoData2$race[which(phenoData2$race == "AA")] <- 1
phenoData2$race[which(phenoData2$race == "Caucasian")] <- 0
phenoData2$race <- as.numeric(as.character(phenoData2$race))

# Optimal matching on a probit PS ----------------------------------------------

m.out1 <- matchit(
  Diabetes ~ age + race + gender + weight + height + Hypertension,
  data = phenoData1,
  method = "optimal",
  distance = "glm",
  link = "probit"
)

m.out2 <- matchit(
  race ~ age + race + gender + weight + height + Diabetes + Hypertension,
  data = phenoData2,
  method = "optimal",
  distance = "glm",
  link = "probit"
)

# Checking balance after full matching -----------------------------------------

summary(m.out1, un = TRUE)
plot(m.out1, type = "jitter", interactive = FALSE)
plot(m.out1,
     type = "density", 
     interactive = FALSE,
     which.xs = ~ age + race + gender + weight + height + Hypertension
)
plot(summary(m.out1))

summary(m.out2, un = TRUE)
plot(m.out2, type = "jitter", interactive = FALSE)
plot(m.out2,
     type = "density", 
     interactive = FALSE,
     which.xs = ~ age  + gender + weight + height + Diabetes + Hypertension
)
plot(summary(m.out2))

# Save matched phenotype data to .csv file -------------------------------------

m.data1 <- match_data(m.out1)
m.data2 <- match_data(m.out2)

# Restore Diabetes status
m.data1$Diabetes[which(m.data1$Diabetes == 1)] <- "Yes"
m.data1$Diabetes[which(m.data1$Diabetes == 0)] <- "No"

# Restore Ethnicity
m.data2$race[which(m.data2$race == 1)] <- "AA"
m.data2$race[which(m.data2$race == 0)] <- "Caucasian"

# Save files
write.csv(m.data1, file = "MAGNet_PhenoData_Matched_Diabetes.csv")
write.csv(m.data2, file = "MAGNet_PhenoData_Matched_Ethnicity.csv")

################################################################################
# CPM CONVERSION                                                               #
################################################################################

# Variable of interest switch --------------------------------------------------

# Change first argument to: 1 for Diabetes, 2 for race
phenoData <- switch(1, m.data1, m.data2)

# Load data --------------------------------------------------------------------

rawCountData <- read.csv("MAGNet_RawCounts.csv", as.is = T, row.names = 1)
gxDataMAGNet <- readRDS("CPMS_SVA_corrected.RDS")
gene_info_all <- readRDS("gene_info_all.RDS") # Ensembl Backup Data file

# Extract protein coding data from Ensembl -------------------------------------

protein_coding <- gene_info_all[gene_info_all$gene_biotype == "protein_coding",]
num_protein_coding_all <- nrow(protein_coding)
#print(num_protein_coding_all)

num_protein_coding <- nrow(protein_coding[rownames(protein_coding) %in% rownames(rawCountData),])
#print(num_protein_coding)

# BLABLABLA --------------------------------------------------------------------

rawCountData <- rawCountData[,colnames(rawCountData) %in% rownames(phenoData)]
cpmData <- cpm(rawCountData, log = FALSE)

keep_genes <- rowSums(cpmData > 1) >= (ncol(cpmData) / 2)
num_genes_kept <- sum(keep_genes)
print(num_genes_kept)

filtered_cpm_matrix <- cpmData[keep_genes, ]

cpmDataMAGNet <- cpmData[rownames(cpmData) %in% rownames(gxDataMAGNet),]

cpmData_log_filtered <- log2(filtered_cpm_matrix + 1)
cpmData_log_magnet <- log2(magnet_cpm + 1)

plot_expression_density <- function(cpmData_log) {
  expression_values <- as.vector(cpmData_log)
  density <- ggplot(data.frame(expression_values), aes(x = expression_values)) +
    geom_density(fill = "skyblue", alpha = 0.6) +
    labs(title = "Expression Value Distribution", x = "Log2(CPM)", y = "Density") +
    theme_minimal() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
  
  plot(density)
}

plot_expression_density(cpmData_log_filtered)
plot_expression_density(cpmData_log_magnet)

sample_tree <- hclust(dist(t(cpmData_log_filtered)), method = "average")
plot(sample_tree, main = "Sample Clustering to Detect Outliers")
# there is one outlier in the data namely P01629, maybe remove? But matching?

# correct for confounding factors, maybe add more? Produces negative values because of shifting but is no problem for WGCNA
design <- model.matrix(~ Library.Pool + gender + Hypertension, data = phenoData1)
adjusted_log_cpm_filtered <- removeBatchEffect(cpmData_log_filtered, covariates = design[, -1])
adjusted_log_cpm_magnet <- removeBatchEffect(cpmData_log_magnet, covariates = design[, -1])

plot_expression_density(adjusted_log_cpm_filtered)
plot_expression_density(adjusted_log_cpm_magnet)

# Save files -------------------------------------------------------------------

saveRDS(adjusted_log_cpm_filtered, "final/gxData_filtered_self.RDS")
saveRDS(adjusted_log_cpm_magnet, "final/gxData_filtered_MAGNet.RDS")

################################################################################