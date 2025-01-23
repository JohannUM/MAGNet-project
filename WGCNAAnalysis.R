# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("RCy3", quietly = TRUE)) BiocManager::install("RCy3")

# Load the necessary libraries ------------------------------------------------

library(here)
library(WGCNA)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RCy3)

allowWGCNAThreads()

# Load the data ------------------------------------------------

setwd(here("data"))

# set the etiology to NF or DCM
etiology <- "DCM"

data_expression <- readRDS(paste(etiology, "/data_", etiology, "_tmm_cpm_log.RDS", sep = ""))
data_samples <- read.csv(paste(etiology, "/data_samples_", etiology, "_Diabetes.csv", sep = ""), row.names = 1)
trait_data <- readRDS(paste(etiology, "/trait_data_", etiology, ".RDS", sep = ""))

# load module_eigengenes, module_labels, module_colors and gene_tree
load(file = paste(etiology, "/network_construction_", etiology, ".RData", sep = ""))

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
names(diabetes) <- "diabetes"
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

# set the module of interest
module <- "tan"
column <- match(module, mod_names)
module_genes <- module_colors == module
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(gene_module_membership[module_genes, column]), abs(gene_trait_significance[module_genes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for Diabetes",
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
)

# get the gene significance and module membership for the module of interest
# look at this to object to see all genes in the module and its significance and membership
gene_significance_membership <- data.frame(
    GeneSignificance = abs(gene_trait_significance[module_genes, 1]),
    ModuleMembership = abs(gene_module_membership[module_genes, column])
)
rownames(gene_significance_membership) <- rownames(gene_module_membership[module_genes, ])

# Perform GO enrichment analysis ------------------------------------------------

genes_module <- rownames(gene_significance_membership)

goResults <- enrichGO(
    gene = genes_module,
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

# Plot TOM heatmap ------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!
# Cant recommend this because it takes a lot of time and need to run WGCNA.R locally to have the TOM object (it is in the gitignore and wont be pushed)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!
# TOM <- readRDS(paste(etiology, "/TOM_", etiology, ".RDS", sep = ""))

# diss_TOM <- 1 - TOM
# plot_TOM <- diss_TOM^7
# diag(plotTOM) <- NA

# sizeGrWindow(15, 15)
# TOMplot(plot_TOM, gene_tree, module_colors, main = "Network heatmap plot, all genes")

# Plot eigengene adjacency ------------------------------------------------

module_eigengenes_trait <- orderMEs(cbind(new_module_eigengenes, diabetes))

# do this in one or next one seperate
sizeGrWindow(11, 11)
par(cex = 0.8)
plotEigengeneNetworks(module_eigengenes_trait, "",
    marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8,
    xLabelsAngle = 90
)


sizeGrWindow(6, 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram",
    marDendro = c(0, 4, 2, 0),
    plotHeatmaps = FALSE
)

par(cex = 1.0)
plotEigengeneNetworks(module_eigengenes_trait, "Eigengene adjacency heatmap",
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = FALSE, xLabelsAngle = 90
)

# Export to Cytoscape ------------------------------------------------

# !!!!!!!!!!! NEED to run WGCNA locally before for this to work !!!!!!!!!!!!!
TOM <- readRDS(paste(etiology, "/TOM_", etiology, ".RDS", sep = ""))

RCy3::cytoscapePing()
if (!grepl("status: Installed", RCy3::getAppStatus("stringApp"))) RCy3::installApp("stringApp")

# Select the module of interest
modules <- c("tan", "darkred")
probes <- names(data_expression)
in_module <- is.finite(match(module_colors, modules))
mod_probes <- probes[in_module]
# Select the corresponding Topological Overlap
mod_TOM <- TOM[in_module, in_module]

dimnames(mod_TOM) <- list(mod_probes, mod_probes)

# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(mod_TOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt", sep = ""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt", sep = ""),
    weighted = TRUE,
    threshold = 0.08,
    nodeNames = mod_probes,
    nodeAttr = module_colors[in_module]
)

edges <- read.table(paste0("CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt"), sep = "")
edges <- edges[, c(1:3)]
colnames(edges) <- c("source", "target", "weight")
RCy3::createNetworkFromDataFrames(edges = edges, title = paste0("module", " ", paste(modules, collapse = "-")))
