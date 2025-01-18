# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")

# Load the necessary libraries ------------------------------------------------

library(here)
library(WGCNA)

allowWGCNAThreads()

# Load the data ------------------------------------------------

setwd(here("data"))

# load the previously filtered and processed data
data_expression <- readRDS("final/data_expression_filtered_magnet.RDS")
data_samples <- read.csv("MAGNet_PhenoData_Matched.csv", row.names = 1)
trait_data <- readRDS("trait_data.RDS")

# load module_eigengenes, module_labels, module_colors and gene_tree
load(file = "network_construction.RData")

data_expression <- t(data_expression) # have to transpose it here to work with the WGCNA functions

# Correlate the module eigengenes with the trait data ------------------------------------------------

n_genes <- ncol(data_expression)
n_samples <- nrow(data_expression)

# recalculate module eigengenes
new_module_eigengenes <- moduleEigengenes(data_expression, module_colors)$eigengenes
module_eigengenes <- orderMEs(new_module_eigengenes)
module_trait_cor <- cor(module_eigengenes, trait_data, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, n_samples)

text_matrix <- paste(signif(module_trait_cor, 2), "\n(",
    signif(module_trait_pvalue, 1), ")",
    sep = ""
)
dim(text_matrix) <- dim(module_trait_cor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
    Matrix = module_trait_cor,
    xLabels = names(trait_data),
    yLabels = names(module_eigengenes),
    ySymbols = names(module_eigengenes),
    colorLabels = FALSE,
    colors = greenWhiteRed(50),
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = paste("Module-trait relationships")
)
