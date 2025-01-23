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
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")

# Load the necessary libraries -------------------------------------------------

library(BiocManager)
library(DESeq2)
library(here)
library(apeglm)
library(Glimma)
library(edgeR)
library(clusterProfiler)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# Load the data ----------------------------------------------------------------

setwd(here("data"))
phenoDataNF <- read.csv("NF/data_samples_NF_Diabetes.csv", row.names = 1)
phenoDataDCM <- read.csv("DCM/data_samples_DCM_Diabetes.csv", row.names = 1)
rawCounts <- read.csv("MAGNet_RawCounts.csv", row.names = 1)

# Extract gene expression data for sample IDs ----------------------------------

ctsNF <- rawCounts[, rownames(phenoDataNF)]
ctsDCM <- rawCounts[, rownames(phenoDataDCM)]

# Construct DESeq DataSet ------------------------------------------------------

ddsNF <- DESeqDataSetFromMatrix(
    countData = ctsNF,
    colData = phenoDataNF,
    design = ~ Library.Pool + RIN + Diabetes
)

ddsDCM <- DESeqDataSetFromMatrix(
    countData = ctsDCM,
    colData = phenoDataDCM,
    design = ~ Library.Pool + RIN + Diabetes
)


# Perform DE Analysis ----------------------------------------------------------

ddsNF <- DESeq(ddsNF, quiet = TRUE)
ddsDCM <- DESeq(ddsDCM, quiet = TRUE)

# Results ----------------------------------------------------------------------

resultsNames(ddsNF)
resultsNames(ddsDCM)

resNF <- results(ddsNF, name = "Diabetes_Yes_vs_No")
resDCM <- results(ddsDCM, name = "Diabetes_Yes_vs_No")

resLFC_NF <- lfcShrink(ddsNF, coef = "Diabetes_Yes_vs_No", type = "apeglm")
resLFC_DCM <- lfcShrink(ddsDCM, coef = "Diabetes_Yes_vs_No", type = "apeglm")

# Plot results -----------------------------------------------------------------

par(mfrow = c(2, 2))

DESeq2::plotMA(resNF, ylim = c(-2, 2))
DESeq2::plotMA(resLFC_NF, ylim = c(-2, 2))

DESeq2::plotMA(resDCM, ylim = c(-2, 2))
DESeq2::plotMA(resLFC_DCM, ylim = c(-2, 2))

EnhancedVolcano(resNF,
    lab = rownames(resNF),
    x = "log2FoldChange",
    y = "pvalue"
)

EnhancedVolcano(resDCM,
    lab = rownames(resDCM),
    x = "log2FoldChange",
    y = "pvalue"
)

# Gene set enrichment analysis NF ----------------------------------------------

resNF <- resNF[order(-resNF$log2FoldChange), ]
gene_list_NF <- resNF$log2FoldChange
names(gene_list_NF) <- rownames(resNF)
gene_list_NF <- gene_list_NF[!is.na(gene_list_NF)]

gsea_go_NF <- gseGO(gene_list_NF,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    eps = 1e-300
)

df_gsea_go_NF <- as.data.frame(gsea_go_NF)

# Gene set enrichment analysis DCM ---------------------------------------------

resDCM <- resDCM[order(-resDCM$log2FoldChange), ]
gene_list_DCM <- resDCM$log2FoldChange
names(gene_list_DCM) <- rownames(resDCM)
gene_list_DCM <- gene_list_DCM[!is.na(gene_list_DCM)]

gsea_go_DCM <- gseGO(gene_list_DCM,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    eps = 1e-300
)

df_gsea_go_DCM <- as.data.frame(gsea_go_DCM)

# Dot Plot ---------------------------------------------------------------------

dotplot(gsea_go_DCM, showCategory = 25, title = "GSE Analysis DCM")
dotplot(gsea_go_NF, showCategory = 20, title = "GSE Analysis NF")

################################################################################
