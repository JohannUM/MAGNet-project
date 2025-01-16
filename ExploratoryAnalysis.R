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

# Load the necessary libraries -------------------------------------------------

library(here)
library(ggplot2)
library(MatchIt)
library(table1)

# Load the data ----------------------------------------------------------------

setwd(here("MAGNet-project/data"))
phenoData <- read.csv("MAGNet_PhenoData_Matched.csv", row.names = 1)
gxData <- readRDS("CPMS_SVA_corrected.RDS")

# Extract gene expression data for matched samples -----------------------------

gxData <- gxData[,rownames(phenoData)]

# Table ------------------------------------------------------------------------

# Set labels for each characteristic and factorize race and diabetes status
label(phenoData$age)     <- "Age"
label(phenoData$height)  <- "Height"
label(phenoData$weight)  <- "Weight"
label(phenoData$race)    <- "Ethnicity"
phenoData$race <- factor(phenoData$race, levels = c("AA","Caucasian"),
                         labels = c("African American","Caucasian"))
phenoData$Diabetes <- factor(phenoData$Diabetes, levels = c(0,1),
                             labels = c("Non-Diabetic","Diabetic"))

# Print the table
table1(~ age + height + weight + race + Hypertension | Diabetes, 
       dat = phenoData,
       footnote = "", 
       overall = TRUE)