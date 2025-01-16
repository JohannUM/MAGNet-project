################################################################################
#                                                                              #
# PRO4002 - Project Period                                                     #
#                                                                              #
# Preprocessing Script                                                         #
#                                                                              #
# K.J. Boessen (i6075128)                                                      #
# N. Imal (i6276013)                                                           #
# J.E. Lottermoser (i6235171)                                                  #
# A. Panghe (i6246854)                                                         #
#                                                                              #
################################################################################

# Install the necessary packages -----------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("MatchIt", quietly = TRUE)) install.packages("MatchIt")


# Load the necessary libraries -------------------------------------------------

library(here)
library(ggplot2)
library(MatchIt)

# Load the phenotype data ------------------------------------------------------

setwd(here("MAGNet-project/data"))
phenoData <- read.csv("MAGNet_PhenoData.csv", row.names = 1)

# Remove Donor, HCM, and PPCM data ---------------------------------------------

phenoData <- subset(phenoData, etiology == "DCM")

# Compute BMI and filter out BMI > 65 ------------------------------------------

phenoData$BMI <- phenoData$weight / ((phenoData$height/100)^2)
phenoData <- subset(phenoData, BMI <= 65)

################################################################################