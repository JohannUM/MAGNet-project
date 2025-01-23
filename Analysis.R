# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

# Load the necessary libraries ------------------------------------------------

library(here)
library(WGCNA)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

allowWGCNAThreads()

# Load the data ------------------------------------------------

setwd(here("data"))

data_expression <- readRDS("DCM/data_DCM_tmm_cpm_log.RDS")
data_samples <- read.csv("DCM/data_samples_DCM_Diabetes.csv", row.names = 1)
trait_data <- readRDS("DCM/trait_data_DCM.RDS")
# Count the number of samples with Diabetes status "Yes"
num_diabetic_yes <- sum(data_samples$Diabetes == "Yes")
print(paste("Number of samples with Diabetes 'Yes':", num_diabetic_yes))
# load module_eigengenes, module_labels, module_colors and gene_tree
load(file = "DCM/network_construction_DCM.RData")

data_expression <- t(data_expression) # have to transpose it here to work with the WGCNA functions

# Correlate the module eigengenes with the trait data ------------------------------------------------

n_genes <- ncol(data_expression)
n_samples <- nrow(data_expression)

# recalculate module eigengenes
new_module_eigengenes <- moduleEigengenes(data_expression, module_colors)$eigengenes
module_eigengenes <- orderMEs(new_module_eigengenes)
module_trait_cor <- cor(module_eigengenes, trait_data, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, n_samples)

sizeGrWindow(9, 5)
text_matrix <- paste(signif(module_trait_cor, 2), "\n(",
    signif(module_trait_pvalue, 1), ")",
    sep = ""
)
dim(text_matrix) <- dim(module_trait_cor)
par(mar = c(6, 8.5, 3, 1))
heatmap <- labeledHeatmap(
    Matrix = module_trait_cor,
    xLabels = names(trait_data),
    yLabels = names(module_eigengenes),
    ySymbols = names(module_eigengenes),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = paste("Module-trait relationships"),
    legendSpace = 0.5
)

# Calculate module membership and gene significance --------------------------------------------------

diabetes <- as.data.frame(trait_data$diabetes)
names(diabates) <- "diabetes"
mod_names <- substring(names(module_eigengenes), 3)
gene_module_membership <- as.data.frame(cor(data_expression, module_eigengenes, use = "p"))
MMP_value <- as.data.frame(corPvalueStudent(as.matrix(gene_module_membership), n_samples))

names(gene_module_membership) <- paste("MM", mod_names, sep = "")
names(MMP_value) <- paste("p.MM", mod_names, sep = "")
gene_trait_significance <- as.data.frame(cor(data_expression, diabetes, use = "p"))
GSP_value <- as.data.frame(corPvalueStudent(as.matrix(gene_trait_significance), n_samples))
names(gene_trait_significance) <- paste("GS.", names(diabetes), sep = "")
names(GSP_value) <- paste("p.GS.", names(diabetes), sep = "")

# Look at module membership vs. gene significance ------------------------------------------------

module <- "lightgreen"
column <- match(module, mod_names)
module_genes <- module_colors == module
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(gene_module_membership[module_genes, column]), abs(gene_trait_significance[module_genes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for Diabetes",
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "darkgreen"
)

genes_paleturquoise <- colnames(data_expression)[module_genes]

goResults <- enrichGO(
    gene = genes_paleturquoise,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

head(goResults)
sizeGrWindow(7, 7)
dotplot(goResults, showCategory = 10, title = "GO Enrichment Analysis")
