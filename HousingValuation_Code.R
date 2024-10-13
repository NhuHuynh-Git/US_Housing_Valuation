---
  title: "US_Housing_Valuation_Predictive_Model"
output: html_notebook
---
  ```{r}
####DATA PREPROCESSING
###DATA TRANSFORMATION
library(readr)
library(dplyr)
HousingValuation <- read.csv("~/Documents/Predictive Analytics/Assessment/HousingValuation.csv")
HousingValuation <- HousingValuation %>% select(-Id)
View(HousingValuation)

##Nominal variables transformation
#LotConfig
lotconfig_Ins <-as.numeric(HousingValuation$LotConfig=="Inside")
lotconfig_Cor <-as.numeric(HousingValuation$LotConfig=="Corner")
lotconfig_Cul <- as.numeric(HousingValuation$LotConfig == "CulDSac")
lotconfig_FR2 <- as.numeric(HousingValuation$LotConfig == "FR2")
lotconfig_FR3 <- as.numeric(HousingValuation$LotConfig == "FR3")

#DwellClass
dwellclass_1Fam <- as.numeric(HousingValuation$DwellClass == "1Fam")
dwellclass_2FmCon <- as.numeric(HousingValuation$DwellClass == "2fmCon")
dwellclass_Duplx <- as.numeric(HousingValuation$DwellClass == "Duplex")
dwellclass_TwnhsE <- as.numeric(HousingValuation$DwellClass == "TwnhsE")
dwellclass_TwnhsI <- as.numeric(HousingValuation$DwellClass == "Twnhs")

#CentralAir
centralAir_Y <- as.numeric(HousingValuation$CentralAir == "Y")
centralAir_N <- as.numeric(HousingValuation$CentralAir == "N")

#GarageType
HousingValuation$GarageType[is.na(HousingValuation$GarageType)] <- "NoGarage"#Replace NA with "NoGarage"

garageType_2Types <- as.numeric(HousingValuation$GarageType == "2Types")
garageType_Attchd <- as.numeric(HousingValuation$GarageType == "Attchd")
garageType_Basment <- as.numeric(HousingValuation$GarageType == "Basment")
garageType_BuiltIn <- as.numeric(HousingValuation$GarageType == "BuiltIn")
garageType_CarPort <- as.numeric(HousingValuation$GarageType == "CarPort")
garageType_Detchd <- as.numeric(HousingValuation$GarageType == "Detchd")
garageType_NoGarage <- as.numeric(HousingValuation$GarageType=="NoGarage")

##Ordinal variables transformation
#LotShape
HousingValuation$LotShape <- factor(HousingValuation$LotShape, levels = 1:4, labels = c(0, 1, 2, 3))
HousingValuation$LotShape <- as.numeric(as.character(HousingValuation$LotShape))

#LandContour
HousingValuation$LandContour <- factor(HousingValuation$LandContour,
                                       levels = c("Lvl","Bnk","HLS","Low"),
                                       labels = c(0, 1, 2, 3))
HousingValuation$LandContour <- as.numeric(as.character(HousingValuation$LandContour))

#Utilities
HousingValuation$Utilities <- factor(HousingValuation$Utilities,
                                     levels = c("AllPub", "NoSewr", "NoSeWa", "ELO"),
                                     labels = c(0, 1, 2, 3))
HousingValuation$Utilities <- as.numeric(as.character(HousingValuation$Utilities))

#Slope
HousingValuation$Slope <- factor(HousingValuation$Slope,
                                 levels = c("Gtl", "Mod", "Sev"),
                                 labels = c(0, 1, 2))
HousingValuation$Slope <- as.numeric(as.character(HousingValuation$Slope))

#OverallQuality
HousingValuation$OverallQuality <- factor(HousingValuation$OverallQuality,
                                          levels = 1:10,
                                          labels = 0:9)
HousingValuation$OverallQuality <- as.numeric(as.character(HousingValuation$OverallQuality))

#OverallCondition
HousingValuation$OverallCondition <- factor(HousingValuation$OverallCondition,
                                            levels = 1:10,
                                            labels = 0:9)
HousingValuation$OverallCondition <- as.numeric(as.character(HousingValuation$OverallCondition))

#ExteriorCondition
HousingValuation$ExteriorCondition <- factor(HousingValuation$ExteriorCondition,
                                             levels = c("Ex", "Gd", "TA", "Fa", "Po"),
                                             labels = c(0, 1, 2, 3, 4))
HousingValuation$ExteriorCondition <- as.numeric(as.character(HousingValuation$ExteriorCondition))

#BasementCondition
HousingValuation$BasementCondition <- factor(HousingValuation$BasementCondition,
                                             levels = c("Ex", "Gd", "TA", "Fa", "Po", "NB"),
                                             labels = c(0, 1, 2, 3, 4, 5))
HousingValuation$BasementCondition <- as.numeric(as.character(HousingValuation$BasementCondition))

#KitchenQuality
HousingValuation$KitchenQuality <- factor(HousingValuation$KitchenQuality,
                                          levels = c("Ex", "Gd", "TA", "Fa", "Po"),
                                          labels = c(0, 1, 2, 3, 4))
HousingValuation$KitchenQuality <- as.numeric(as.character(HousingValuation$KitchenQuality))

#PavedDrive
HousingValuation$PavedDrive <- factor(HousingValuation$PavedDrive,
                                      levels = c("Y", "P", "N"),
                                      labels = c(0, 1, 2))
HousingValuation$PavedDrive <- as.numeric(as.character(HousingValuation$PavedDrive))

#Combine nominal and ordinal values, and remove original columns dataset
HousingValuation <- cbind(HousingValuation,lotconfig_Ins,lotconfig_Cor, lotconfig_Cul,lotconfig_FR2, lotconfig_FR3, dwellclass_1Fam, dwellclass_2FmCon,dwellclass_Duplx, dwellclass_TwnhsE, dwellclass_TwnhsI, centralAir_Y, centralAir_N, garageType_2Types, garageType_Attchd, garageType_Basment,garageType_BuiltIn,garageType_CarPort, garageType_Detchd, garageType_NoGarage)

library(dplyr)
HousingValuation <- HousingValuation %>%
  select(-LotConfig, -DwellClass, -CentralAir, -GarageType)

View(HousingValuation) 
summary(HousingValuation
)
```

```{r}
###SUMMARY STATISTICS
##Numeric variables
#LotArea
cat("LotArea:\n")
cat("Mean: ", mean(HousingValuation$LotArea, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$LotArea, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$LotArea, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$LotArea, na.rm = TRUE), "\n\n")

#TotalBSF
cat("TotalBSF:\n")
cat("Mean: ", mean(HousingValuation$TotalBSF, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$TotalBSF, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$TotalBSF, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$TotalBSF, na.rm = TRUE), "\n\n")

#LowQualFinSF
cat("LowQualFinSF:\n")
cat("Mean: ", mean(HousingValuation$LowQualFinSF, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$LowQualFinSF, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$LowQualFinSF, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$LowQualFinSF, na.rm = TRUE), "\n\n")

#LivingArea
cat("LivingArea:\n")
cat("Mean: ", mean(HousingValuation$LivingArea, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$LivingArea, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$LivingArea, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$LivingArea, na.rm = TRUE), "\n\n")

#OpenPorchSF
cat("OpenPorchSF:\n")
cat("Mean: ", mean(HousingValuation$OpenPorchSF, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$OpenPorchSF, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$OpenPorchSF, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$OpenPorchSF, na.rm = TRUE), "\n\n")

#PoolArea
cat("PoolArea:\n")
cat("Mean: ", mean(HousingValuation$PoolArea, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$PoolArea, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$PoolArea, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$PoolArea, na.rm = TRUE), "\n\n")

#SalePrice
cat("SalePrice:\n")
cat("Mean: ", mean(HousingValuation$SalePrice, na.rm = TRUE), "\n")
cat("Median: ", median(HousingValuation$SalePrice, na.rm = TRUE), "\n")
cat("Max: ", max(HousingValuation$SalePrice, na.rm = TRUE), "\n")
cat("Standard Deviation: ", sd(HousingValuation$SalePrice, na.rm = TRUE), "\n\n")

##Nonimal variables:
cat("Nominal Variables Counts:\n")
cat("LotConfig - Inside: ", sum(lotconfig_Ins, na.rm = TRUE), "\n")
cat("LotConfig - Corner: ", sum(lotconfig_Cor, na.rm = TRUE), "\n")
cat("LotConfig - CulDSac: ", sum(lotconfig_Cul, na.rm = TRUE), "\n")
cat("LotConfig - FR2: ", sum(lotconfig_FR2, na.rm = TRUE), "\n")
cat("LotConfig - FR3: ", sum(lotconfig_FR3, na.rm = TRUE), "\n")
cat("DwellClass - 1Fam: ", sum(dwellclass_1Fam, na.rm = TRUE), "\n")
cat("DwellClass - 2FmCon: ", sum(dwellclass_2FmCon, na.rm = TRUE), "\n")
cat("DwellClass - Duplx: ", sum(dwellclass_Duplx, na.rm = TRUE), "\n")
cat("DwellClass - TwnhsE: ", sum(dwellclass_TwnhsE, na.rm = TRUE), "\n")
cat("DwellClass - TwnhsI: ", sum(dwellclass_TwnhsI, na.rm = TRUE), "\n")
cat("CentralAir - Y: ", sum(centralAir_Y, na.rm = TRUE), "\n")
cat("CentralAir - N: ", sum(centralAir_N, na.rm = TRUE), "\n")
cat("GarageType - 2Types: ", sum(garageType_2Types, na.rm = TRUE), "\n")
cat("GarageType - Attchd: ", sum(garageType_Attchd, na.rm = TRUE), "\n")
cat("GarageType - Basment: ", sum(garageType_Basment, na.rm = TRUE), "\n")
cat("GarageType - BuiltIn: ", sum(garageType_BuiltIn, na.rm = TRUE), "\n")
cat("GarageType - CarPort: ", sum(garageType_CarPort, na.rm = TRUE), "\n")
cat("GarageType - Detchd: ", sum(garageType_Detchd, na.rm = TRUE), "\n")
cat("GarageType - NoGarage: ", sum(garageType_NoGarage, na.rm = TRUE), "\n\n")

##Ordinal variables:
cat("Ordinal Variables Counts:\n")
# LotShape
cat("LotShape:\n")
lotshape_counts <- table(HousingValuation$LotShape)
cat("Reg (0): ", lotshape_counts["0"], "\n")
cat("IR1 (1): ", lotshape_counts["1"], "\n")
cat("IR2 (2): ", lotshape_counts["2"], "\n")
cat("IR3 (3): ", lotshape_counts["3"], "\n\n")

# LandContour
cat("LandContour:\n")
landcontour_counts <- table(HousingValuation$LandContour)
cat("Lvl (0): ", landcontour_counts["0"], "\n")
cat("Bnk (1): ", landcontour_counts["1"], "\n")
cat("HLS (2): ", landcontour_counts["2"], "\n")
cat("Low (3): ", landcontour_counts["3"], "\n\n")

# Utilities
cat("Utilities:\n")
utilities_counts <- table(HousingValuation$Utilities)
cat("AllPub (0): ", utilities_counts["0"], "\n")
cat("NoSewr (1): ", utilities_counts["1"], "\n")
cat("NoSeWa (2): ", utilities_counts["2"], "\n")
cat("ELO (3): ", utilities_counts["3"], "\n\n")

# Slope
cat("Slope:\n")
slope_counts <- table(HousingValuation$Slope)
cat("Gtl (0): ", slope_counts["0"], "\n")
cat("Mod (1): ", slope_counts["1"], "\n")
cat("Sev (2): ", slope_counts["2"], "\n\n")

# OverallQuality
cat("OverallQuality:\n")
overallquality_counts <- table(HousingValuation$OverallQuality)
for (i in 0:9) {
  cat(i, ": ", overallquality_counts[as.character(i)], "\n")
}
cat("\n")

# OverallCondition
cat("OverallCondition:\n")
overallcondition_counts <- table(HousingValuation$OverallCondition)
for (i in 0:9) {
  cat(i, ": ", overallcondition_counts[as.character(i)], "\n")
}
cat("\n")

# ExteriorCondition
cat("ExteriorCondition:\n")
exteriorcondition_counts <- table(HousingValuation$ExteriorCondition)
cat("Ex (0): ", exteriorcondition_counts["0"], "\n")
cat("Gd (1): ", exteriorcondition_counts["1"], "\n")
cat("TA (2): ", exteriorcondition_counts["2"], "\n")
cat("Fa (3): ", exteriorcondition_counts["3"], "\n")
cat("Po (4): ", exteriorcondition_counts["4"], "\n\n")

# BasementCondition
cat("BasementCondition:\n")
basementcondition_counts <- table(HousingValuation$BasementCondition)
cat("Ex (0): ", basementcondition_counts["0"], "\n")
cat("Gd (1): ", basementcondition_counts["1"], "\n")
cat("TA (2): ", basementcondition_counts["2"], "\n")
cat("Fa (3): ", basementcondition_counts["3"], "\n")
cat("Po (4): ", basementcondition_counts["4"], "\n")
cat("NB (5): ", basementcondition_counts["5"], "\n\n")

# KitchenQuality
cat("KitchenQuality:\n")
kitchenquality_counts <- table(HousingValuation$KitchenQuality)
cat("Ex (0): ", kitchenquality_counts["0"], "\n")
cat("Gd (1): ", kitchenquality_counts["1"], "\n")
cat("TA (2): ", kitchenquality_counts["2"], "\n")
cat("Fa (3): ", kitchenquality_counts["3"], "\n")
cat("Po (4): ", kitchenquality_counts["4"], "\n\n")

# PavedDrive
cat("PavedDrive:\n")
paveddrive_counts <- table(HousingValuation$PavedDrive)
cat("Y (0): ", paveddrive_counts["0"], "\n")
cat("P (1): ", paveddrive_counts["1"], "\n")
cat("N (2): ", paveddrive_counts["2"], "\n")
```

```{r}
###VARIABLE DISTRIBUTION INVESTIGATION
##Histograms and summary statistics of continuous variables
#Histograms
library(purrr) 
library(tidyr) 
library(ggplot2)
library(dplyr)
continuous <- HousingValuation %>% 
  select(LotArea, TotalBSF, LowQualFinSF, LivingArea, OpenPorchSF, PoolArea, SalePrice)

continuous %>%
  keep(is.numeric) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") + 
  geom_histogram() +
  labs(title = "Histograms of Continuous Variables",
       x = "Value",
       y = "Frequency")

#Summary statistics
summary(continuous)

```
```{r}
###HANDLING MISSING VALUES
summary(HousingValuation)

cat("Summary statistics of LivingArea_Original: \n")
summary(HousingValuation$LivingArea)
cat("Summary statistics of YearBuilt_Original: \n")
summary(HousingValuation$YearBuilt)
library(tidyverse)

##Method 1: fill missing values with median for YearBuilt and mean for LivingArea
HousingValuation_mean <- HousingValuation

HousingValuation_mean$YearBuilt[is.na(HousingValuation_mean$YearBuilt)] <-median(HousingValuation$YearBuilt, na.rm = TRUE)
HousingValuation_mean$LivingArea[is.na(HousingValuation_mean$LivingArea)] <-mean(HousingValuation$LivingArea, na.rm = TRUE)
HousingValuation_mean

cat("Summary statistics of LivingArea_Method 1: \n")
summary(HousingValuation_mean$LivingArea)
cat("Summary statistics of YearBuilt_Method 1: \n")
summary(HousingValuation_mean$YearBuilt)

#LivingArea
plot(density(HousingValuation_mean$LivingArea), col = "red", 
     main = "LivingArea: Original (Blue) vs Mean Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$LivingArea, na.rm = TRUE), col = "blue")

#YearBuilt
plot(density(HousingValuation_mean$YearBuilt), col = "red", 
     main = "YearBuilt: Original (Blue) vs Mean Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$YearBuilt, na.rm = TRUE), col = "blue")


##Method 2: deleting all missing values
HousingValuation_delete <- HousingValuation[complete.cases(HousingValuation),] 

#LivingArea
plot(density(HousingValuation_delete$LivingArea), col = "red", 
     main = "LivingArea: Original (Blue) vs Deletion Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$LivingArea, na.rm = TRUE), col = "blue")

#YearBuilt
plot(density(HousingValuation_delete$YearBuilt), col = "red", 
     main = "YearBuilt: Original (Blue) vs Deletion Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$YearBuilt, na.rm = TRUE), col = "blue")

cat("Summary statistics of LivingArea_Method 2: \n")
summary(HousingValuation_delete$LivingArea)
cat("Summary statistics of YearBuilt_Method 2: \n")
summary(HousingValuation_delete$YearBuilt)


##Method 3: replacing missing values with 0
HousingValuation_zeros <- HousingValuation
HousingValuation_zeros[is.na(HousingValuation_zeros)] <- 0

plot(density(HousingValuation_zeros$LivingArea), col = "red", 
     main = "LivingArea: Original (Blue) vs Zero Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$LivingArea, na.rm = TRUE), col = "blue")

#YearBuilt
plot(density(HousingValuation_zeros$YearBuilt), col = "red", 
     main = "YearBuilt: Original (Blue) vs Zero Transformation (Red)",
     xlab = "Value", ylab = "Density")
lines(density(HousingValuation$YearBuilt, na.rm = TRUE), col = "blue")

cat("Summary statistics of LivingArea_Method 3: \n")
summary(HousingValuation_zeros$LivingArea)
cat("Summary statistics of YearBuilt_Method 3: \n")
summary(HousingValuation_zeros$YearBuilt)
```

```{r}
###MULTICOLLINEARITY CHECK AND DIMENSION REDUCTION
##Call all relevant packages
library(Amelia) 
library(psych)
library(corrplot)
library(caret) 
library(plotly) 
library(ggplot2) 
library(GGally) 

cor.plot(HousingValuation_mean, numbers = TRUE)

##Dimension reduction
#Remove target variable
target <- HousingValuation_mean$SalePrice
HousingValuation_mean_subset <- subset(HousingValuation_mean, select=-c(SalePrice))

ggcorr(HousingValuation_mean_subset, label=TRUE)

#Create data matrix
matrix <- data.matrix(HousingValuation_mean_subset)
corr_HousingValuation_mean <- cor(matrix)

#Find the highly correlated variables with a cutoff of 0.5
highlyCorr_HousingValuation_mean <- findCorrelation(corr_HousingValuation_mean, cutoff=0.5)
names(HousingValuation_mean_subset)[highlyCorr_HousingValuation_mean]

#Carry out dimension reduction_1
subset_selected_1 <- subset(HousingValuation_mean_subset, select = -c(OverallQuality, FullBath, TotalRmsAbvGrd))

ggcorr(subset_selected_1, label=TRUE)

#Merge target back
subset_selected_1$SalePrice <-target
View(subset_selected_1)

##Distribution of selected variables
library(purrr) 
library(tidyr) 
library(ggplot2) 
library(dplyr)
subset_selected_1 %>%
  gather() %>%  
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") + 
  geom_histogram() 


#Replace 0s in the dataset with 0.00001
subset_selected_1 <- subset_selected_1 %>%
  mutate(
    LotArea = ifelse(LotArea == 0, 0.00001, LotArea),
    OpenPorchSF = ifelse(OpenPorchSF == 0, 0.00001, OpenPorchSF),
    TotalBSF = ifelse(TotalBSF == 0, 0.00001, TotalBSF)
  )

#Transform variables with high skewness
trans_LotArea <- log(subset_selected_1$LotArea)
trans_OpenPorchSF <- log(subset_selected_1$OpenPorchSF)
trans_TotalBSF<- log(subset_selected_1$TotalBSF)

#Plot histograms for original and transformed distributions
par(mfrow = c(3, 2))

#LotArea
hist(subset_selected_1$LotArea, col = "orange", main = "Original LotArea")
hist(trans_LotArea, col = "orange", main = "Transformed LotArea")

#OpenPorchSF
hist(subset_selected_1$OpenPorchSF, col = "orange", main = "Original OpenPorchSF")
hist(trans_OpenPorchSF, col = "orange", main = "Transformed OpenPorchSF")

#TotalBSF
hist(subset_selected_1$TotalBSF, col = "orange", main = "Original TotalBSF")
hist(trans_TotalBSF, col = "orange", main = "Transformed TotalBSF")

cols <-c('LotArea', 'OpenPorchSF', 'TotalBSF')
subset_selected_1[(cols)] <- log(subset_selected_1[(cols)])


```

```{r}
####PREDICTIVE MODEL BUIDLING
###REGRESSION MODELS

# Set the sample size (2/3 of the data)
smp_size <- floor(2/3 * nrow(subset_selected_1))  
set.seed(2) 

# Shuffle the dataset randomly
subset_selected_1 <- 
  subset_selected_1[sample(nrow(subset_selected_1)), ] 

# Split the data into training and testing sets
subset_selected_1.train <- subset_selected_1[1:smp_size, ]   
subset_selected_1.test <- subset_selected_1[(smp_size+1):nrow(subset_selected_1), ]   

# Define the formula and fit the model
formula = SalePrice ~ .
model <- lm(formula = formula, data = subset_selected_1.train)
summary(model)$coefficients

# Generate the equation from the coefficients
as.formula( 
  paste0("y ~ ", round(coefficients(model)[1], 2), " + ",  
         paste(sprintf("%.2f * %s", coefficients(model)[-1], 
                       names(coefficients(model)[-1])),  
               collapse=" + ")
  )
)

# Predict SalePrice for the test set
subset_selected_1.test$predicted.SalePrice <- predict(model, newdata = subset_selected_1.test)

# Display the actual and predicted values
print("Actual Values:") 
head(subset_selected_1.test$SalePrice, 5)

print("Predicted Values:") 
head(subset_selected_1.test$predicted.SalePrice, 5)

# Plot actual vs predicted values
pl1 <- subset_selected_1.test %>%  
  ggplot(aes(SalePrice, predicted.SalePrice)) + 
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(colour = 'red')) + 
  xlab('Actual value of SalePrice') + 
  ylab('Predicted value of SalePrice') + 
  theme_bw() 
ggplotly(pl1)

# Calculate the error (residuals)
error <- subset_selected_1.test$SalePrice - subset_selected_1.test$predicted.SalePrice

# Calculate RMSE
rmse <- sqrt(mean(error^2))
print(paste("Root Mean Square Error: ", rmse))
```

```{r}
##Regression model 1
subset_selected_2 <- subset_selected_1 %>% select(-LotArea,-OpenPorchSF,-TotalBSF)

# Set the sample size (2/3 of the data)
smp_size <- floor(2/3 * nrow(subset_selected_2))  
set.seed(2) 

# Shuffle the dataset randomly
subset_selected_2 <- 
  subset_selected_2[sample(nrow(subset_selected_2)), ] 

# Split the data into training and testing sets
subset_selected_2.train <- subset_selected_2[1:smp_size, ]   
subset_selected_2.test <- subset_selected_2[(smp_size+1):nrow(subset_selected_2), ]   

# Define the formula and fit the model
formula = SalePrice ~ .
model <- lm(formula = formula, data = subset_selected_2.train)
summary(model)$coefficients

# Generate the equation from the coefficients
as.formula( 
  paste0("y ~ ", round(coefficients(model)[1], 2), " + ",  
         paste(sprintf("%.2f * %s", coefficients(model)[-1], 
                       names(coefficients(model)[-1])),  
               collapse=" + ")
  )
)

# Predict SalePrice for the test set
subset_selected_2.test$predicted.SalePrice <- predict(model, newdata = subset_selected_2.test)

# Display the actual and predicted values
print("Actual Values:") 
head(subset_selected_2.test$SalePrice, 5)

print("Predicted Values:") 
head(subset_selected_2.test$predicted.SalePrice, 5)

# Plot actual vs predicted values
pl1 <- subset_selected_2.test %>%  
  ggplot(aes(SalePrice, predicted.SalePrice)) + 
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(colour = 'red')) + 
  xlab('Actual value of SalePrice') + 
  ylab('Predicted value of SalePrice') + 
  theme_bw() 
ggplotly(pl1)

# Calculate the error (residuals)
error <- subset_selected_2.test$SalePrice - subset_selected_2.test$predicted.SalePrice

# Calculate RMSE
rmse <- sqrt(mean(error^2))
print(paste("Root Mean Square Error: ", rmse))

```
```{r}
##Regression model 2
subset_selected_3 <- subset_selected_1 %>% select(-OpenPorchSF,-TotalBSF)

# Set the sample size (2/3 of the data)
smp_size <- floor(2/3 * nrow(subset_selected_3))  
set.seed(2) 

# Shuffle the dataset randomly
subset_selected_3 <- 
  subset_selected_3[sample(nrow(subset_selected_3)), ] 

# Split the data into training and testing sets
subset_selected_3.train <- subset_selected_3[1:smp_size, ]   
subset_selected_3.test <- subset_selected_3[(smp_size+1):nrow(subset_selected_3), ]   

# Define the formula and fit the model
formula = SalePrice ~ .
model <- lm(formula = formula, data = subset_selected_3.train)
summary(model)$coefficients

# Generate the equation from the coefficients
as.formula( 
  paste0("y ~ ", round(coefficients(model)[1], 2), " + ",  
         paste(sprintf("%.2f * %s", coefficients(model)[-1], 
                       names(coefficients(model)[-1])),  
               collapse=" + ")
  )
)

# Predict SalePrice for the test set
subset_selected_3.test$predicted.SalePrice <- predict(model, newdata = subset_selected_3.test)

# Display the actual and predicted values
print("Actual Values:") 
head(subset_selected_3.test$SalePrice, 5)

print("Predicted Values:") 
head(subset_selected_3.test$predicted.SalePrice, 5)

# Plot actual vs predicted values
pl1 <- subset_selected_3.test %>%  
  ggplot(aes(SalePrice, predicted.SalePrice)) + 
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(colour = 'red')) + 
  xlab('Actual value of SalePrice') + 
  ylab('Predicted value of SalePrice') + 
  theme_bw() 
ggplotly(pl1)

# Calculate the error (residuals)
error <- subset_selected_3.test$SalePrice - subset_selected_3.test$predicted.SalePrice

# Calculate RMSE
rmse <- sqrt(mean(error^2))
print(paste("Root Mean Square Error: ", rmse))

```

```{r}
##Regression model 3
subset_selected_4 <- subset_selected_1 %>% select(-LivingArea,-YearBuilt)

# Set the sample size (2/3 of the data)
smp_size <- floor(2/3 * nrow(subset_selected_4))  
set.seed(2) 

# Shuffle the dataset randomly
subset_selected_4 <- 
  subset_selected_4[sample(nrow(subset_selected_4)), ] 

# Split the data into training and testing sets
subset_selected_4.train <- subset_selected_4[1:smp_size, ]   
subset_selected_4.test <- subset_selected_4[(smp_size+1):nrow(subset_selected_4), ]   

# Define the formula and fit the model
formula = SalePrice ~ .
model <- lm(formula = formula, data = subset_selected_4.train)
summary(model)$coefficients

# Generate the equation from the coefficients
as.formula( 
  paste0("y ~ ", round(coefficients(model)[1], 2), " + ",  
         paste(sprintf("%.2f * %s", coefficients(model)[-1], 
                       names(coefficients(model)[-1])),  
               collapse=" + ")
  )
)

# Predict SalePrice for the test set
subset_selected_4.test$predicted.SalePrice <- predict(model, newdata = subset_selected_4.test)

# Display the actual and predicted values
print("Actual Values:") 
head(subset_selected_4.test$SalePrice, 5)

print("Predicted Values:") 
head(subset_selected_4.test$predicted.SalePrice, 5)

# Plot actual vs predicted values
pl1 <- subset_selected_4.test %>%  
  ggplot(aes(SalePrice, predicted.SalePrice)) + 
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(colour = 'red')) + 
  xlab('Actual value of SalePrice') + 
  ylab('Predicted value of SalePrice') + 
  theme_bw() 
ggplotly(pl1)

# Calculate the error (residuals)
error <- subset_selected_4.test$SalePrice - subset_selected_4.test$predicted.SalePrice

# Calculate RMSE
rmse <- sqrt(mean(error^2))
print(paste("Root Mean Square Error: ", rmse))

```

```{r}
###DECISION TREES
#Install and load necessary libraries
#install.packages('rpart.plot', dependencies = TRUE) 
#install.packages('rpart', dependencies = TRUE) 
library(corrplot) 
library(GGally) 
library(rpart) 
library(rpart.plot)

#Remove target variable
target <- HousingValuation_mean$SalePrice
HousingValuation_mean_tree <- subset(HousingValuation_mean, select=-c(SalePrice))

#Create data matrix
matrix <- data.matrix(HousingValuation_mean_tree)
corr_HousingValuation_mean_tree <- cor(matrix)

#Find the highly correlated variables with a cutoff of 0.5
highlyCorr_HousingValuation_mean_tree <- findCorrelation(corr_HousingValuation_mean_tree, cutoff=0.5)
names(HousingValuation_mean_tree)[highlyCorr_HousingValuation_mean_tree]

#Carry out dimension reduction
#Remove FullBath and TotalRmsAbvGrd
tree_selected <- subset(HousingValuation_mean_tree, select = -c(FullBath, TotalRmsAbvGrd))

#Merge target back
tree_selected$SalePrice <- target
View(tree_selected)

#Split the data into training and test sets
smp_size <- floor(2/3 * nrow(tree_selected))  
set.seed(2) 
tree.selected <- tree_selected[sample(nrow(tree_selected)), ] 
tree_selected.train <- tree.selected[1:smp_size, ]   
tree_selected.test <- tree.selected[(smp_size+1):nrow(tree.selected), ]  

##Decision tree model 1
#Define the formula for decision tree
formula <- SalePrice ~ . 

#Build the decision tree using rpart
dtree <- rpart(formula, data=tree_selected.train, method="anova")
print(dtree$variable.importance)
rpart.plot(dtree, type = 4, fallen.leaves = FALSE)

#Print the decision tree
print(dtree)

#Make predictions on the test set
predicted.SalePrice <- predict(dtree, tree_selected.test)

#Compare actual vs predicted values
print("Actual Values") 
head(tree_selected.test$SalePrice[1:5]) 
print("Predicted Values") 
head(predicted.SalePrice[1:5])

#Calculate RMSE for the unpruned tree
error <- tree_selected.test$SalePrice - predicted.SalePrice
rmse <- sqrt(mean(error^2)) 
print(paste("Root Mean Square Error (Unpruned): ", rmse))

#Print complexity parameter table and find best cp value
printcp(dtree)
best_cp <- dtree$cptable[which.min(dtree$cptable[,"xerror"]),"CP"] 
print(paste("Best CP: ", best_cp))

##Decision tree model 2
#Prune the tree with cp = 0.011628
pruned_dtree_1 <- prune(dtree, cp = 0.011628)  
rpart.plot(pruned_dtree_1, type = 4, fallen.leaves = FALSE)

#Make predictions with the pruned tree
predicted_pruned_1.SalePrice <- predict(pruned_dtree_1, tree_selected.test)
error_pruned_1 <- tree_selected.test$SalePrice - predicted_pruned_1.SalePrice
rmse_pruned_1 <- sqrt(mean(error_pruned_1^2))  # Corrected RMSE calculation
print(paste("Root Mean Square Error (Pruned 0.011628): ", rmse_pruned_1))

##Decision tree model 3
#Prune the tree with cp = 0.02
pruned_dtree_2 <- prune(dtree, cp = 0.02)
rpart.plot(pruned_dtree_2, type = 4, fallen.leaves = FALSE)

#Make predictions with the second pruned tree
predicted_pruned_2.SalePrice <- predict(pruned_dtree_2, tree_selected.test)
error_pruned_2 <- tree_selected.test$SalePrice - predicted_pruned_2.SalePrice
rmse_pruned_2 <- sqrt(mean(error_pruned_2^2))  # Corrected RMSE calculation
print(paste("Root Mean Square Error (Pruned 0.02): ", rmse_pruned_2))
```

