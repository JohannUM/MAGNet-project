################################################################################
#                                                                              #
# PRO4002 - Project Period                                                     #
#                                                                              #
# Sample Matching Script                                                       #
#                                                                              #
# K.J. Boessen (i6075128)                                                      #
# N. Imal (i6276013)                                                           #
# J.E. Lottermoser (i6235171)                                                  #
# A. Panghe (i6246854)                                                         #
#                                                                              #
################################################################################

# Install the necessary packages -----------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("MatchIt", quietly = TRUE)) install.packages("MatchIt")
if (!requireNamespace("rlemon", quietly = TRUE)) install.packages("rlemon")
if (!requireNamespace("optmatch", quietly = TRUE)) install.packages("optmatch")


# Load the necessary libraries -------------------------------------------------

library(here)
library(MatchIt)
library(rlemon)
library(optmatch)

# Load the phenotype data ------------------------------------------------------

setwd(here("MAGNet-project/data"))
phenoData <- read.csv("MAGNet_PhenoData_Processed.csv", row.names = 1)

# Convert Diabetes to numeric --------------------------------------------------

phenoData$Diabetes[which(phenoData$Diabetes == "Yes")] <- 1
phenoData$Diabetes[which(phenoData$Diabetes == "No")] <- 0
phenoData$Diabetes = as.numeric(as.character(phenoData$Diabetes)) 

# Optimal matching on a probit PS ----------------------------------------------

m.out <- matchit(Diabetes ~ age + race + gender + weight + height + 
                    Hypertension,
                  data = phenoData,
                  method = "optimal",
                  distance = "glm",
                  link = "probit")

# Checking balance after full matching -----------------------------------------

summary(m.out, un = FALSE)
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out, type = "density", interactive = FALSE,
     which.xs = ~age + race + gender + weight + height + Hypertension)
plot(summary(m.out))

# Save matched phenotype data to .csv file -------------------------------------

m.data <- match_data(m.out)
write.csv(m.data, file = "MAGNet_PhenoData_Matched.csv")

################################################################################