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
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("Glimma", quietly = TRUE)) BiocManager::install("Glimma")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")

# Load the necessary libraries -------------------------------------------------

library(BiocManager)
library(DESeq2)
library(here)
library(apeglm)
library(Glimma)
library(edgeR)
library(clusterProfiler)

# Load the data ----------------------------------------------------------------

setwd(here("data"))
phenoData0 <- read.csv("MAGNet_PhenoData_DCM.csv", row.names = 1)
phenoData1 <- read.csv("MAGNet_PhenoData_Matched_Diabetes.csv", row.names = 1)
phenoData2 <- read.csv("MAGNet_PhenoData_Matched_Ethnicity.csv", row.names = 1)
rawCounts <- read.csv("MAGNet_RawCounts.csv", row.names = 1)

# Convert ethnicity to factors -------------------------------------------------

phenoData2$race <- as.factor(phenoData2$race)

# Extract gene expression data for matched samples -----------------------------

cts0 <- rawCounts[,rownames(phenoData0)]
cts1 <- rawCounts[,rownames(phenoData1)]

# Construct DESeq DataSet ------------------------------------------------------

dds0 <- DESeqDataSetFromMatrix(countData = cts0,
                               colData = phenoData0,
                               design = ~ Library.Pool + Diabetes)

dds1 <- DESeqDataSetFromMatrix(countData = cts1,
                               colData = phenoData1,
                               design = ~ Library.Pool + Diabetes)


# Perform DE Analysis ----------------------------------------------------------

dds0 <- DESeq(dds0, quiet = TRUE)
dds1 <- DESeq(dds1, quiet = TRUE)

# Results ----------------------------------------------------------------------

resultsNames(dds0)
resultsNames(dds1)

res0 <- results(dds0, name = "Diabetes_Yes_vs_No")
res1 <- results(dds1, name = "Diabetes_Yes_vs_No")

resLFC0 <- lfcShrink(dds0, coef = "Diabetes_Yes_vs_No", type = "apeglm")
resLFC1 <- lfcShrink(dds1, coef = "Diabetes_Yes_vs_No", type = "apeglm")

# Plot results -----------------------------------------------------------------

par(mfrow=c(2,2))

DESeq2::plotMA(res0, ylim = c(-2,2))
DESeq2::plotMA(resLFC0, ylim = c(-2,2))

DESeq2::plotMA(res1, ylim = c(-2,2))
DESeq2::plotMA(resLFC1, ylim = c(-2,2))

# Gene set enrichment analysis ALL ---------------------------------------------

logFC0 <- res0$log2FoldChange
names(logFC0) <- rownames(res0)
logFC0 <- sort(logFC0, decreasing = TRUE)

# Up-regulated genes

upreg0 <- rownames(res0)[res0$pvalue < 0.05 & res0$log2FoldChange > 0]

# GO gene set enrichment analysis
gsea_go0 <- enrichGO(upreg0, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENSEMBL", 
                     ont = "BP", 
                     universe = rownames(res0))

df_gsea_go0 <- as.data.frame(gsea_go0)

# Gene set enrichment analysis Diabetes ----------------------------------------

logFC1 <- res1$log2FoldChange
names(logFC1) <- rownames(res1)
logFC1 <- sort(logFC1, decreasing = TRUE)

# Up-regulated genes
upreg1 <- rownames(res1)[res1$pvalue < 0.05 & res1$log2FoldChange > 0]

# GO gene set enrichment analysis
gsea_go1 <- enrichGO(upreg1, 
                 OrgDb = org.Hs.eg.db, 
                 keyType = "ENSEMBL", 
                 ont = "BP", 
                 universe = rownames(res1))

df_gsea_go1 <- as.data.frame(gsea_go1)

################################################################################