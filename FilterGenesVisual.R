# Install the necessary packages and set options ------------------------------------------------
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")

# Load the necessary libraries ------------------------------------------------

library(here)
library(ggplot2)

# Load the data ------------------------------------------------

setwd(here("data"))

# Take preprocessed data here
data_expression <- readRDS("CPMS_SVA_corrected.RDS")
data_samples <- read.csv("phenoData.csv", as.is = T, row.names = 1)
data_exon_lengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, row.names = 1)

if (!all(rownames(data_expression) == rownames(data_exon_lengths))) {
  print("Row names of data_expression and data_exon_lengths do not match or are not in the same order.")
}

data_samples <- data_samples[data_samples$etiology == "DCM", ]
data_expression <- data_expression[, rownames(data_samples)]

# log transform the data
data_log <- log2(data_expression + 1)

expression_values <- as.vector(data_log) # Flatten the matrix to a vector
density <- ggplot(data.frame(expression_values), aes(x = expression_values)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  labs(title = "Expression Value Distribution", x = "Log2(CPM)", y = "Density") +
  theme_minimal() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))

plot(density)

# Filter rows where the row mean is above 2
filtered_data_log <- data_log[rowMeans(data_log, na.rm = TRUE) > 2, ]

# Print the dimensions of the filtered data
print(dim(filtered_data_log))

saveRDS(filtered_data_log, "data_expression_filtered.RDS")
saveRDS(data_samples, "data_samples_filtered.RDS")
