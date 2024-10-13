# Predictive Models on Property Prices

## Table of Contents
-	[Project Overview](#project-overview)
-	[Data Sources](#data-sources)
-	[Tools](#tools)
-	[Data Exploration](#data-exploration)
-	[Data Preprocessing](#data-preprocessing)
-	[Predictive Model Building](#predictive-model-building)
-	[Results and Findings](#results-and-findings)
-	[Business Recommendations](#business-recommendations)
-	[Limitations](#limitations)
-	[References](#references)

### Project Overview
This project focuses on developing predictive models for a property assessor’s office in the USA to enhance housing price assessments. In particular, the project encompasses data preprocessing, model building and evaluation to select the most effective predictive model for real-world application.

### Data Sources
Housing Price Data: The main dataset for training the predictive model in this project is the “HousingValuation.csv” file, which includes over 3,900 property sales records with 31 features for each property.

### Tools
-	Excel – Data Exploration
-	R – Data Preprocessing, Predictive Model Building 

### Data Exploration
Before preprocessing and analysis, the following exploration tasks are performed:
1.	External research on domain understanding
2.	Data loading and inspection
3.	Metadata and variable understanding

### Data Preprocessing
In the phase of preprocessing, the following tasks are performed:
1.	Variable (nominal and numeric) transformation 
2.	Variable distribution investigation
3.	Handling missing values 
4.	Multicollinearity check and dimension reduction

### Predictive Model Building
-	Regression modelling (3 models with feature selection)
-	Decision tree modelling (3 models with pruning)

### Results and findings
Using RMSE as an accuracy metric, the results of model comparison can be summarized as follows
-	Regression model 2, after removing OpenPorchSF and TotalBSF for feature selection purpose, is more optimal that model 1 and 3
-	Decision tree 1, with OverallQuality as the first split, is shown to perform better than model 2 and 3 with lower RMSE
-	Comparing these two optimal models, in this particular dataset, using the regression model is better in terms of producing the outcome with the smaller errors

### Business recommendations:
Based on the analysis, the recommendations are as follows
-	The decision tree excels in interpretability, segmenting properties based on key attributes like OverallQuality, LivingArea, and TotalBSF (Total Basement Square Feet). This makes it easy for stakeholders to understand how different property characteristics impact pricing. 
-	On the other hand, the regression model (after removing OpenPorchSF and TotalBSF) provides a continuous equation for predicting property prices, offering precise, numeric insights into the impact of each variable. 
-	Given that RMSE was used as the performance metric, and the regression model showed slightly lower errors, it may be the better choice for price prediction accuracy in this case. 

### Limitations:
-	Time-Bound Dataset:
The dataset covers sales from 2006 to 2010, and the housing market dynamics might have changed since then, limiting the applicability of the models for future predictions without retraining on updated data.
-	Impact of Feature Selection:
Removing features like OpenPorchSF and TotalBSF from the regression model might have affected predictive performance. Important variables may have been excluded, leading to suboptimal results.
-	Limited Scope of Metrics (RMSE Only):
RMSE, while useful for measuring prediction error, does not provide insights into other aspects such as model interpretability, robustness, or bias-variance trade-off. Using additional metrics like MAE or R² could provide a more comprehensive evaluation.

### References:
1. Gomez, J. (2019, March 27). 8 critical factors that influence a home’s value | Opendoor. Opendoor. https://www.opendoor.com/w/blog/factors-that-influence-home-value
2. Gordon Linoff and Michael Berry, Data Mining Techniques, 3rd edition, Wiley, 2011
3. Martin, E. J. (n.d.). How much is my house worth? Bankrate. https://www.bankrate.com/real-estate/how-much-is-my-house-worth/

