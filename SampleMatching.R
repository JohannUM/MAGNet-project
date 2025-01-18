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


# Load the necessary libraries -------------------------------------------------

library(here)
library(MatchIt)

# Load the phenotype data ------------------------------------------------------

setwd(here("data"))
phenoData <- read.csv("MAGNet_PhenoData_Processed.csv", row.names = 1)
phenoData1 <- phenoData
phenoData2 <- phenoData

# Convert Diabetes to numeric --------------------------------------------------

<<<<<<< HEAD
phenoData1$Diabetes[which(phenoData1$Diabetes == "Yes")] <- 1
phenoData1$Diabetes[which(phenoData1$Diabetes == "No")] <- 0
phenoData1$Diabetes = as.numeric(as.character(phenoData1$Diabetes)) 

# Optimal matching on a probit PS ----------------------------------------------

m.out1 <- matchit(Diabetes ~ age + race + gender + weight + height + 
                    Hypertension,
                  data = phenoData1,
                  method = "optimal",
                  distance = "glm",
                  link = "probit")

# Checking balance after full matching -----------------------------------------

summary(m.out1, un = FALSE)
plot(m.out1, type = "jitter", interactive = FALSE)
plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~age + race + gender + weight + height + Hypertension)
plot(summary(m.out1))

# Save matched phenotype data to .csv file -------------------------------------

m.data1 <- match_data(m.out1)
write.csv(m.data1, file = "MAGNet_PhenoData_Matched_Diabetes.csv")

################################################################################

# Convert Diabetes to numeric --------------------------------------------------

phenoData2$race[which(phenoData2$race == "AA")] <- 1
phenoData2$race[which(phenoData2$race == "Caucasian")] <- 0
phenoData2$race = as.numeric(as.character(phenoData2$race)) 

# Optimal matching on a probit PS ----------------------------------------------

m.out2 <- matchit(race ~ age + gender + weight + height + Diabetes + 
                   Hypertension,
                 data = phenoData2,
                 method = "optimal",
                 distance = "glm",
                 link = "probit")

# Checking balance after full matching -----------------------------------------

summary(m.out2, un = FALSE)
plot(m.out2, type = "jitter", interactive = FALSE)
plot(m.out2, type = "density", interactive = FALSE,
     which.xs = ~age + gender + weight + height + Diabetes + Hypertension)
plot(summary(m.out2))

# Save matched phenotype data to .csv file -------------------------------------

m.data2 <- match_data(m.out2)
write.csv(m.data2, file = "MAGNet_PhenoData_Matched_Ethnicity.csv")
=======
phenoData$Diabetes[which(phenoData$Diabetes == "Yes")] <- 1
phenoData$Diabetes[which(phenoData$Diabetes == "No")] <- 0
phenoData$Diabetes <- as.numeric(as.character(phenoData$Diabetes))

# Optimal matching on a probit PS ----------------------------------------------

m.out <- matchit(
     Diabetes ~ age + race + gender + weight + height +
          Hypertension,
     data = phenoData,
     method = "optimal",
     distance = "glm",
     link = "probit"
)

# Checking balance after full matching -----------------------------------------

summary(m.out, un = TRUE)
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out,
     type = "density", interactive = FALSE,
     which.xs = ~ age + race + gender + weight + height + Hypertension
)
plot(summary(m.out))

# Save matched phenotype data to .csv file -------------------------------------

m.data <- match_data(m.out)
m.data$Diabetes[which(m.data$Diabetes == 1)] <- "Yes"
m.data$Diabetes[which(m.data$Diabetes == 0)] <- "No"

write.csv(m.data, file = "MAGNet_PhenoData_Matched.csv")
>>>>>>> Cleanup

################################################################################
