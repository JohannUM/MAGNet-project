################################################################################
#                                                                              #
# PRO4002 - Project Period                                                     #
#                                                                              #
# Exploratory Analysis Script                                                  #
#                                                                              #
# K.J. Boessen (i6075128)                                                      #
# N. Imal (i6276013)                                                           #
# J.E. Lottermoser (i6235171)                                                  #
# A. Panghe (i6246854)                                                         #
#                                                                              #
################################################################################

# Install the necessary packages -----------------------------------------------

if (!requireNamespace("table1", quietly = TRUE)) install.packages("table1")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidy")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pcaMethods", quietly = TRUE)) install.packages("pcaMethods")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("rapportools", quietly = TRUE)) install.packages("rapportools")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")


# Load the necessary libraries -------------------------------------------------

library(here)
library(ggplot2)
library(MatchIt)
library(table1)
library(pcaMethods)


# Load the data ----------------------------------------------------------------

setwd(here("data"))
phenoData0 <- read.csv("DCM/data_samples_DCM_Diabetes.csv", row.names = 1)
phenoData1 <- read.csv("MAGNet_PhenoData_Matched_Diabetes.csv", row.names = 1)
phenoData2 <- read.csv("MAGNet_PhenoData_Matched_Ethnicity.csv", row.names = 1)
gxData <- readRDS("DCM/data_DCM_tmm_cpm_log.RDS")

# Extract gene expression data for matched samples -----------------------------

gxData0 <- gxData[,rownames(phenoData0)]
gxData1 <- gxData[,rownames(phenoData1)]
gxData2 <- gxData[,rownames(phenoData2)]

# Function ---------------------------------------------------------------------

exploreData <- function (phenoData, gxData, strat) {
  
  # Set labels for each characteristic and factorize race and diabetes status
  label(phenoData$age)     <- "Age"
  label(phenoData$height)  <- "Height"
  label(phenoData$weight)  <- "Weight"
  label(phenoData$race)    <- "Ethnicity"
  phenoData$race <- factor(phenoData$race, levels = c("AA","Caucasian"),
                            labels = c("African American","Caucasian"))
  phenoData$Diabetes <- factor(phenoData$Diabetes, levels = c("No","Yes"),
                                labels = c("NonDiabetic","Diabetic"))
  
  # Print a table
  table <- switch(strat, 
         "diabetes" = table1(~ age + height + weight + race + Hypertension | Diabetes, dat = phenoData),
         "race" = table1(~ age + height + weight + Diabetes + Hypertension | race, dat = phenoData),
         "DCM" = table1(~ age + height + weight + Diabetes + race + Hypertension | etiology, dat = phenoData)
         )
  table
  
  # Perform PCA on the gene expression data for 10 principal components
  pcaRes <- pca(t(gxData), nPcs = 10)

  # Plot PCA figure showing gender, age, and diabetes 
  plotData <- cbind(data.frame(pcaRes@scores), phenoData)
  
  plot <- switch(strat, 
         "diabetes" = ggplot(plotData, aes(x = PC1, y = PC2)) + geom_point(aes(color = Diabetes, shape = gender, size = age)),
         "race" = ggplot(plotData, aes(x = PC1, y = PC2)) + geom_point(aes(color = race, shape = gender, size = age)),
         "DCM" = ggplot(plotData, aes(x = PC1, y = PC2)) + geom_point(aes(color = Diabetes, shape = gender, size = age))
         )
  plot
  
}

# Test function

exploreData(phenoData0,gxData0,"DCM")
exploreData(phenoData1,gxData1,"diabetes")
exploreData(phenoData2,gxData2,"race")


pcaRes <- pca(t(gxData0), nPcs = 10)
plotData <- cbind(data.frame(pcaRes@scores), phenoData0)
plot <- ggplot(plotData, aes(x = PC1, y = PC2)) + geom_point(aes(color = Diabetes, shape = gender, size = age))
plot
plot <- ggplot(plotData, aes(x = PC1, y = PC2)) + geom_point(aes(color = Library.Pool, shape = Diabetes))
plot

# Numeric Plots ----------------------------------------------------------------

numeric_data <- phenoData0[, sapply(phenoData0, is.numeric)]

# Check structure of numeric data
str(numeric_data)

# Boxplot for all numeric columns
boxplot(numeric_data, main = "Boxplot for Numeric Variables", las = 2)