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
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

# Load the necessary libraries -------------------------------------------------

library(BiocManager)
library(DESeq2)
library(here)

# Load the data ----------------------------------------------------------------

setwd(here("data"))
phenoData <- read.csv("MAGNet_PhenoData_Matched.csv", row.names = 1)
rawCounts <- read.csv("MAGNet_RawCounts.csv", row.names = 1)

# Extract gene expression data for matched samples -----------------------------

cts <- rawCounts[,rownames(phenoData)]

# Construct DESeq DataSet ------------------------------------------------------

# Factorize diabetes status
phenoData$Diabetes <- factor(phenoData$Diabetes, levels = c(0,1),
                             labels = c("Non-Diabetic","Diabetic"))

# Create DataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = phenoData,
                              design= ~ Library.Pool + Diabetes)

# Perform DE Analysis ----------------------------------------------------------

dds <- DESeq(dds)

# Results

resultsNames(dds) # lists the coefficients
res <- results(dds, name="Diabetes_Diabetic_vs_Non.Diabetic")
