# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages()("here")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr")

# Load the necessary libraries ------------------------------------------------

library(here)
library(WGCNA)
library(dplyr)

allowWGCNAThreads()

# Load the data ------------------------------------------------

setwd(here("data"))

# load the previously filtered and processed data
data_expression <- readRDS("NF/data_NF_tmm_cpm_log.RDS")
data_samples <- read.csv("NF/data_samples_NF_Matched_Diabetes.csv", row.names = 1)

gsg <- goodSamplesGenes(data_expression, verbose = 3)
print(gsg$allOK)

data_expression <- t(data_expression) # have to transpose it here to work with the WGCNA functions

# Create trait data transformation ------------------------------------------------
trait_data <- data.frame(
    sex = ifelse(is.na(data_samples$gender), NA, ifelse(data_samples$gender == "Female", 0, 1)),
    ethnicity = ifelse(is.na(data_samples$race), NA, ifelse(data_samples$race == "AA", 1, 0)),
    age = data_samples$age,
    weight = data_samples$weight,
    height = data_samples$height,
    diabetes = ifelse(is.na(data_samples$Diabetes), NA, ifelse(data_samples$Diabetes == "Yes", 1, 0)),
    hypertension = ifelse(is.na(data_samples$Hypertension), NA, ifelse(data_samples$Hypertension == "Yes", 1, 0))
)
trait_data <- mutate_all(trait_data, function(x) as.numeric(as.character(x)))
rownames(trait_data) <- rownames(data_samples)

collectGarbage()

if (!all(rownames(trait_data) == rownames(data_expression))) {
    print("The rownames of trait_data and data_expression do not match.")
}
saveRDS(trait_data, "trait_data.RDS")

sample_tree <- hclust(dist(data_expression), method = "average")
trait_colors <- numbers2colors(trait_data, signed = FALSE)

plotDendroAndColors(sample_tree, trait_colors,
    groupLabels = names(trait_data), cex.dendroLabels = 0.5,
    main = "Sample dendrogram and trait heatmap DCM"
)

legend("topright",
    legend = c("Female", "Caucasian"),
    fill = c("white", "white"),
    border = c("black", "black", "black", "black"), cex = 0.8
)

# Pick soft threshold ------------------------------------------------
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

sft <- pickSoftThreshold(data_expression, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# picked the soft threshold 5 because of the scale free topology fit
soft_power <- 5

# Network construction and module detection ------------------------------------------------
adjacency <- adjacency(data_expression, power = soft_power)

TOM <- TOMsimilarity(adjacency)
diss_TOM <- 1 - TOM

saveRDS(TOM, "TOM_NF.RDS")

gene_tree <- hclust(as.dist(diss_TOM), method = "average")

plot(gene_tree,
    xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04
)

min_module_size <- 30

dynamic_modules <- cutreeDynamic(
    dendro = gene_tree, distM = diss_TOM,
    deepSplit = 2, pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
)
table(dynamic_modules)


dynamic_colors <- labels2colors(dynamic_modules)

plotDendroAndColors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors"
)

# Merging modules ------------------------------------------------

ME_list <- moduleEigengenes(data_expression, colors = dynamic_colors)
module_eigengenes <- ME_list$eigengenes

ME_dissimilarity <- 1 - cor(module_eigengenes)

ME_tree <- hclust(as.dist(ME_dissimilarity), method = "average")

plot(ME_tree,
    main = "Clustering of module eigengenes",
    xlab = "", sub = ""
)

ME_diss_thres <- 0.25

abline(h = ME_diss_thres, col = "red")

merge <- mergeCloseModules(data_expression, dynamic_colors, cutHeight = ME_diss_thres, verbose = 3)

merged_colors <- merge$colors

merged_ME <- merge$newMEs

plotDendroAndColors(gene_tree, cbind(dynamic_colors, merged_colors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
)

colors <- c("grey", standardColors(60))
module_colors <- merged_colors
module_eigengenes <- merged_ME
module_labels <- match(module_colors, colors) - 1

save(module_eigengenes, module_labels, module_colors, gene_tree, file = "NF/network_construction_NFv2.RData")
