# Practical assignment 3
# Dr Reza Rafiee, 29th Oct. 2019
# Analysis: multiple imputation modelling using Amelia and mice packages in R
# The input csv file includes numeric values with some missing data
# The aim is (but not limited) to the following items:
#   1) Read the CSV file including 220 samples with 17 features, 
#   2) Assess the fraction of missing across all features and samples,
#   3) Use multiple imputation techniques to impute the missingness (both Amelia and mice packages),
#   4) Visualise the imputation results (including imputed values vs. observed, etc.) ,
#   5) Finally comparing the obtained results from two techniques using a scatter plot.
#-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Libraries will be in this section
#
library(Amelia)
library(mice)
#----------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     User functions will be in this section 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Chagning the range of input matrix into [0 1]
Data_Range_Into_01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################################################################
#                                   Main programme - Start
#####################################################################################################
# Set working directory and reading the CSV file including 220 samples with 17 features.
getwd()
setwd("D:/Live")  # change the path to your working directory including the following csv file 
# Read my input csv file from the working directory
Input_dataset <- as.data.frame(t(read.csv("PA3_CSC3062_With_Missing_RR_2019_220.csv",row.names = 1)))
class(Input_dataset)  # [1] "data.frame"
attributes(Input_dataset)
dim(Input_dataset) # [1] 220  17
#----------------------------------------------------------------------------------------------------
#                     First analysis: assess the fraction of missing across all features and samples
#----------------------------------------------------------------------------------------------------
# Running multiple imputation using Amelia package (m=20)
Number_of_imputation <- 20
datasets_models_using_amelia <- amelia(x = Input_dataset, m = Number_of_imputation, p2s = 1, frontend = FALSE) 
summary(datasets_models_using_amelia)
# Fraction Missing
# feature_1       0.009090909
# feature_2       0.027272727
# feature_3       0.000000000
# feature_4       0.018181818
# feature_5       0.013636364
# feature_6       0.013636364
# feature_7       0.018181818
# feature_8       0.027272727
# feature_9       0.031818182
# feature_10      0.013636364
# feature_11      0.022727273
# feature_12      0.018181818
# feature_13      0.027272727
# feature_14      0.013636364
# feature_15      0.009090909
# feature_16      0.009090909
# feature_17      0.009090909

######################################################################################## 
# Adding element wise of the complete datasets obtained from amelia
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
Sum_matrix_2 <- t(datasets_models_using_amelia$imputations[[1]])  # The first complete dataset
for (i_index in 2:Number_of_imputation) {
  Sum_matrix_2 <- Sum_matrix_2 + t(datasets_models_using_amelia$imputations[[i_index]])
}  # for
# Find the average
avg_matrix_elementwis_2 <- Sum_matrix_2/Number_of_imputation

Final_imputed_dataset_using_amelia <- as.data.frame(t(avg_matrix_elementwis_2))
  
min(avg_matrix_elementwis_2)  # what is the min of element in the final matrix: -0.21617
max(avg_matrix_elementwis_2)   # what is the max of element in the final matrix: 1.047615

#Final_imputed_dataset_using_amelia <- Data_Range_Into_01(Final_imputed_dataset_using_amelia) 
######################################################################################## 
# Imputation using mice package
imputated_by_mice <- mice(Input_dataset,m=Number_of_imputation, maxit=Number_of_imputation,seed = 12345) 
# Now aggregating the complete datasets obtained from mice
Sum_matrix_3 <- complete(imputated_by_mice,1)   # The first imputed dataset obtained from mice package
# Adding 
for (i in 2:Number_of_imputation) { 
  Sum_matrix_3 <- Sum_matrix_3 + complete(imputated_by_mice,i)
}
Sum_matrix_3 <- Sum_matrix_3/Number_of_imputation
Final_imputed_dataset_using_mice <- as.data.frame(Sum_matrix_3)

# Till here, we have imputed two datasets obtained from Amelia and mice packages
# Now, try to compare the results 
matplot(Final_imputed_dataset_using_mice, Final_imputed_dataset_using_amelia, type="p")
######################################################################################## 

####################################################################################################
#                                   Main programme - END
#####################################################################################################
