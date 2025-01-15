install.packages("ggplot2")
library(ggplot2)
# Function to find outliers
find_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  outliers <- which(x < lower_bound | x > upper_bound)
  return(outliers)
}
# Read csv file
data <- read.csv("C:/Users/Alexandra/Desktop/MSB Project sem 1/MAGNET_SampleData_18112022.csv")
# print(data)
# Summary statistics
summary(data)
# Display column names
# colnames(data)

# # Inspect the first few rows
# head(data)
# 
# # Check the structure of the data
# str(data)
# 
# # Replace 'ColumnName' with the actual column name
# column_data <- data$ColumnName
# 
# # Check for missing values
# sum(is.na(column_data))

# Convert to numeric if it's not already
numeric_column <- as.numeric(column_data)
# Remove NA values
cleaned_data <- na.omit(numeric_column)

# Check the cleaned data
summary(cleaned_data)

# matrix_data <- as.matrix(data)
# print(matrix_data)
# Filter numeric columns
numeric_data <- data[, sapply(data, is.numeric)]

# Check structure of numeric data
str(numeric_data)

# Boxplot for all numeric columns
boxplot(numeric_data, main = "Boxplot for Numeric Columns", las = 2)

# Find outliers
# Apply the function to each numeric column
outliers_list <- lapply(numeric_data, find_outliers)

# Print outliers
outliers_list

# Histogram for better distribution understanding
hist(numeric_data$age, main = "Histogram of Age", xlab = "Age", breaks = 20)

# weight and height
# Print the height column (optional)
print(data$height)

# Convert height from cm to meters (if needed)
data$height_m <- data$height / 100

# Compute BMI
data$BMI <- data$weight / (data$height_m^2)

# Filter out BMI > 75
data <- data[data$BMI <= 75, ]

# Inspect the updated data
summary(data$BMI)
head(data)

# Check for rows with missing BMI
missing_bmi <- data[is.na(data$BMI), ]

# Print rows with missing values
print(missing_bmi)

# Summary statistics for BMI
summary(data$BMI)

# Histogram to visualize BMI distribution
hist(data$BMI, 
     main = "Histogram of BMI", 
     xlab = "BMI", 
     col = "lightblue", 
     breaks = 20)

#Diabetes
# Create a table of counts for Diabetes
diabetes_counts <- table(data$Diabetes)

# Bar plot
barplot(diabetes_counts, main = "Bar Plot of Diabetes", xlab = "Diabetes", ylab = "Count", col = "pink")

# Unique values in Diabetes
unique(data$Diabetes)

# Filter rows with unexpected values
outliers_diabetes <- data[!data$Diabetes %in% c("No", "Yes"), ]
print(outliers_diabetes)

# Extract the count of people with "Yes" for Diabetes
count_with_diabetes <- diabetes_counts["Yes"]
print(paste("Number of people with diabetes:", count_with_diabetes))

# remove instances of nf (controls), ethyology nf and hcm, maybe pbcm
subset_data <- data[!(data$etiology %in% c("nf", "hcm", "pbcm")), ]
print(subset_data)
summary(subset_data)
# Create a table of counts for Diabetes
diabetes_counts <- table(data$Diabetes)

# Bar plot
barplot(diabetes_counts, main = "Bar Plot of Diabetes", xlab = "Diabetes", ylab = "Count", col = "pink")

# Unique values in Diabetes
unique(data$Diabetes)

# Filter rows with unexpected values
outliers_diabetes <- data[!data$Diabetes %in% c("No", "Yes"), ]
print(outliers_diabetes)

# Extract the count of people with "Yes" for Diabetes
count_with_diabetes <- diabetes_counts["Yes"]
print(paste("Number of people with diabetes:", count_with_diabetes))
count_with_diabetes <- diabetes_counts["No"]
print(paste("Number of people with diabetes:", count_with_diabetes))
