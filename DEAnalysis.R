################################################################################
#                                                                              #
# PRO4002 - Project Period                                                     #
#                                                                              #
# Differential Expression Analysis Script                                      #
#                                                                              #
# K.J. Boessen (i6075128)                                                      #
# N. Imal (i6276013)                                                           #
# J.E. Lottermoser (i6235171)                                                  #
# A. Panghe (i6246854)                                                         #
#                                                                              #
################################################################################

# Install the necessary packages -----------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")


# Load the necessary libraries -------------------------------------------------

library(BiocManager)
library(DESeq2)
library(here)
library(apeglm)

# Load the data ----------------------------------------------------------------

setwd(here("data"))
phenoData1 <- read.csv("MAGNet_PhenoData_Matched_Diabetes.csv", row.names = 1)
phenoData2 <- read.csv("MAGNet_PhenoData_Matched_Ethnicity.csv", row.names = 1)
rawCounts <- read.csv("MAGNet_RawCounts.csv", row.names = 1)

# Convert ethnicity to factors -------------------------------------------------

phenoData2$race <- as.factor(phenoData2$race)

# Extract gene expression data for matched samples -----------------------------

cts1 <- rawCounts[,rownames(phenoData1)]
cts2 <- rawCounts[,rownames(phenoData2)]

# Construct DESeq DataSet ------------------------------------------------------

dds1 <- DESeqDataSetFromMatrix(countData = cts1,
                               colData = phenoData1,
                               design = ~ Library.Pool + Diabetes)
dds2 <- DESeqDataSetFromMatrix(countData = cts2,
                               colData = phenoData2,
                               design = ~ Library.Pool + race)

# Perform DE Analysis ----------------------------------------------------------

dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# Results ----------------------------------------------------------------------

resultsNames(dds1)
resultsNames(dds2)

res1 <- results(dds1, name = "Diabetes_Yes_vs_No")
res2 <- results(dds2, name = "race_Caucasian_vs_AA")

resLFC1 <- lfcShrink(dds1, coef = "Diabetes_Yes_vs_No", type = "apeglm")
resLFC2 <- lfcShrink(dds2, coef = "race_Caucasian_vs_AA", type = "apeglm")

# Plot results -----------------------------------------------------------------

par(mfrow=c(2,2))

DESeq2::plotMA(res1, ylim = c(-2,2))
DESeq2::plotMA(resLFC1, ylim = c(-2,2))

DESeq2::plotMA(res2, ylim = c(-2,2))
DESeq2::plotMA(resLFC2, ylim = c(-2,2))

#plotCounts(dds1, gene = which.min(res1$padj), intgroup = "Diabetes")
#plotCounts(dds2, gene = which.min(res1$padj), intgroup = "race")

################################################################################
