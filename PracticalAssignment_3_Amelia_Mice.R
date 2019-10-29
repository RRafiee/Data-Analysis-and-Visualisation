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
#                     Addition of the two input matrices element-wise  
add_element_wise <- function(x1_matrix,y1_matrix) {
                          Sum_matrix <- x1_matrix   # the sum of two matrix will be here
                          for (i in 1:nrow(x1_matrix))  {  # i is the index of row
                            for (j in 1:ncol(x1_matrix)) {
                              Sum_matrix[i,j] <- Sum_matrix[i,j] + y1_matrix[i,j] 
                            }  # for j
                          }  # for i
                          return(Sum_matrix)  # return the sum of two matrices 
} # function add_element_wise
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Chagning the range of input matrix into [0 1]
Data_Range_Into_01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Imputation using mice package 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MultipleImputationModeling <- function(Input_Sample)
{
  # Columns must be features and rows must be samples
  # Show the missing data pattern
  md.pattern(Input_Sample)
  # Column should be variables (features or gene)
  # maxit: a scalar giving the number of iterations. The default is 5
  # m: Number of multiple imputations. 
  NumofIteration <- 5 
  imputated_by_mice <- mice(Input_Sample,m=20,maxit=NumofIteration,seed = 123)  
  #summary(imputated_by_mice)
  ### blue is observed, red is imputed
  #stripplot(imputated_by_mice, PLG_S5, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #stripplot(imputated_by_mice, MMP3_S2, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #stripplot(imputated_by_mice, RGS8_S5, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #imputated_by_mice2 <- mice.mids(imputated_by_mice)
  # verification
  #identical(imputated_by_mice$imp, imputated_by_mice2$imp)
  ### density plot of head circumference per imputation
  ### blue is observed, red is imputed
  #densityplot(imputated_by_mice, ~PLG_S5|.imp)
  #densityplot(imputated_by_mice, ~MMP3_S2|.imp)
  #densityplot(imputated_by_mice, ~RGS8_S5|.imp)
  #xyplot(imputated_by_mice,...)
  #==========================================
  X <- complete(imputated_by_mice, action = "long", include = TRUE)[, -2]
  test <- as.mids(X, where = NULL, .imp = ".imp", .id = ".id")
  is.mids(test)
  Test_dat <- complete(test, action = "long", include = TRUE)
  #================================================
  #output1 <- matrix(unlist(imputated_by_mice2), ncol = 17, byrow = TRUE)
  #imputated_by_mice$data[1:17]
  #firstSample <- imputated_by_mice$data[3,1:17]
  ColumnNoPlus2MoreVaraibles <- ncol(Input_Sample) + 2
  Imputed_mat <- matrix(nrow = nrow(Input_Sample) ,ncol = ColumnNoPlus2MoreVaraibles ,0.0) # ncol = 19
  nrow(Test_dat)
  #for (impx in 0:5) # number of iterations
  for(idx in 1:nrow(Input_Sample)) # number of samples
  {
    for (j in 3:ColumnNoPlus2MoreVaraibles)   # number of columns
    {
      if (is.na(Test_dat[idx,j]))
      {
        sum.imp <- 0
        for (impx in 1:NumofIteration+1) # number of iterations: 6
        {
          sum.imp <- Test_dat[nrow(Input_Sample)*impx+idx,j] + sum.imp
        }
        Imputed_mat[idx,j] <- sum.imp/(NumofIteration+1)
      }
      else
      {
        Imputed_mat[idx,j] <- Test_dat[idx,j]
      }
    }
  }
  Imputed_mat <- Imputed_mat[,-c(1:2)]
  Imputed_mat <- t(Imputed_mat)
  colnames(Imputed_mat) <- rownames(Input_Sample)
  rownames(Imputed_mat) <- colnames(Input_Sample)
  return(Imputed_mat)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################################################################
#                                   Main programme - Start
#####################################################################################################
# Set working directory and reading the CSV file including 220 samples with 17 features.
getwd()
setwd("D:/Live")
# Read my input csv file from working directory
Input_dataset <- as.data.frame(t(read.csv("PA3_CSC3062_With_Missing_RR_2019_220.csv",row.names = 1)))
class(Input_dataset)  # [1] "data.frame"
attributes(Input_dataset)
dim(Input_dataset) # [1] 220  17
#----------------------------------------------------------------------------------------------------
#                     First analysis: assess the fraction of missing across all features and samples
#----------------------------------------------------------------------------------------------------
# Running multiple imputation using Amelia package (m=20)
Number_of_imputation <- 5
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

Sum_matrix <- t(datasets_models_using_amelia$imputations[[1]])  # The first complete dataset
for (i_index in 2:Number_of_imputation) {
  Sum_matrix <- add_element_wise(Sum_matrix,t(datasets_models_using_amelia$imputations[[i_index]]))
} # for
avg_matrix_elementwise <- Sum_matrix/Number_of_imputation

min(avg_matrix_elementwise) # what is the min of element in the final matrix: -0.1179939
max(avg_matrix_elementwise) # what is the max of element in the final matrix: 1.061036

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# or using the following for loop without a function

Sum_matrix_2 <- t(datasets_models_using_amelia$imputations[[1]])  # The first complete dataset
for (i_index in 2:Number_of_imputation) {
  Sum_matrix_2 <- Sum_matrix_2 + t(datasets_models_using_amelia$imputations[[i_index]])
}  # for
# Find the average
avg_matrix_elementwis_2 <- Sum_matrix_2/Number_of_imputation

Final_imputed_dataset_using_amelia <- as.data.frame(t(avg_matrix_elementwis_2))
  
min(avg_matrix_elementwis_2)  # what is the min of element in the final matrix: -0.1179939
max(avg_matrix_elementwise)   # what is the max of element in the final matrix: 1.061036
######################################################################################## 
# Call the function for multiple imputation using mice 
Final_imputed_dataset_using_mice <- MultipleImputationModeling(Input_dataset)  # Numeric dataset
Final_imputed_dataset_using_mice <- as.data.frame(t(Final_imputed_dataset_using_mice)) # as dataframe
# Till here, we have imputed two datasets obtained from Amelia and mice
# Now, try to compare the results 


######################################################################################## 
####################################################################################################
#                                   Main programme - END
#####################################################################################################
