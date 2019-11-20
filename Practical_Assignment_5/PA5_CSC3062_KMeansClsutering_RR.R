# Practical assignment 5: consensus clustering using k_means algorithm
# Clustering analysis
# Dr Reza Rafiee, 20th Nov. 2019
# The input csv file includes numeric values without missingness
# The aim is (but not limited) to the following items:
#   1) Read an input CSV file including 220 samples and 17 features from input.
#   This is a complete dataset without any missing that you can download from the link below
#   https://github.com/RRafiee/Data-Analysis-and-Visualisation.
#   2) Run k-means clustering algorithm on this dataset and identify the cluster labels (single run).
#   3) Try to use k-means clustering in multiple runs with different initialisation settings and 
#   evaluate and compare the results of multiple runs of the clustering algorithm
#   4) Combine the cluster labels (optional).
#   5) Visualise the outcome of the clustering results (from single run or aggregated)
#   using PCA when samples are coloured based on your clustering labels .
#   Use the following colours in link with a sample subgroup.
#   1: blue; 2: red; 3: yellow; 4:green
# 
# All right reserved!
#-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Libraries will be in this section
library(FactoMineR)
library(factoextra)
library(pca3d)
library(parallel) # For mclapply() that speeds up analysis using multiple cpu cores
library(cluster)  # daisy() for Dissimilarity Matrix Calculation
#----------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     User functions will be in this section 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Changing the range of input matrix into [0 1]
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
# Read an input csv file from the working directory
Complete_dataaset_220 <- as.data.frame(read.csv("Complete_Dataset_220_17_CSC3062_RR_2019.csv",row.names = 1))
min(Complete_dataaset_220, na.rm = TRUE) # 0
max(Complete_dataaset_220, na.rm = TRUE) # 0.9994831
class(Complete_dataaset_220)  # [1] "data.frame"
attributes(Complete_dataaset_220)
dim(Complete_dataaset_220) # [1]  17 220
#----------------------------------------------------------------------------------------------------
# 2) Run k-means clustering algorithm on this dataset and identify the cluster labels (single run)
#----------------------------------------------------------------------------------------------------
K_means_Model <- kmeans(t(Complete_dataaset_220),centers = 4, iter.max = 50,nstart = 10) # nstart: how many random sets should be chosen?
plot(t(Complete_dataaset_220),col=K_means_Model$cluster) # 2D and only two features
points(K_means_Model$centers, col = 1:4, pch = 10, cex = 5)  # cluster centers
# Visualising the cluster labels of clustering results 
# using a PCA method,fviz_cluster() {factoextra library}
fviz_cluster(K_means_Model,t(Complete_dataaset_220),geom = "point",shape = NULL,labelsize = 8,ellipse.type = "norm")

#----------------------------------------------------------------------------------------------------
# 3) Try to use k-means clustering in multiple runs with different initialisation settings
#----------------------------------------------------------------------------------------------------
# Creating a matrix of all cluster labels when running k-means multiple times  
set.seed(20)
Number_of_runs <- 4
Matrix_labels_different_runs <- matrix(nrow = ncol(Complete_dataaset_220), ncol = Number_of_runs,0)
rownames(Matrix_labels_different_runs) <- colnames(Complete_dataaset_220)
# --------- Using a for loop ----------
for (j in 1:Number_of_runs) {
  #j <- 1
  pdf(paste("KmeansClusterPlot_",j,".pdf",sep = "")) 
  K_means_Model <- kmeans(t(Complete_dataaset_220),centers = 4, iter.max = 50,nstart = 5) #
  #Using several random starts (nstart> 1) is often recommended.
  Matrix_labels_different_runs[,j] <- K_means_Model$cluster
  # Visualising the cluster labels of clustering results 
  # using a PCA method,fviz_cluster() {factoextra library}
  fviz_cluster(K_means_Model,t(Complete_dataaset_220),geom = "point",shape = NULL,labelsize = 8,ellipse.type = "norm")
} # for
dev.off()
# --------- Using lapply() ----------
# Using lapply() function instead of the above for loop
K_means_lapply <- lapply(1:Number_of_runs,
                         function(i) 
                           kmeans(t(Complete_dataaset_220),centers = 4, iter.max = 50,nstart = 5))

for (j in 1:Number_of_runs) {
  j <- 2
  Matrix_labels_different_runs [,j]<- K_means_lapply[[j]]$cluster
  fviz_cluster(K_means_lapply[[j]],t(Complete_dataaset_220),
               geom = "point",shape = NULL,labelsize = 8,ellipse.type = "norm")
} # for 

# --------- Using mclapply() when using Linux or Mac OS ----------
K_means_mclapply <- mclapply(1:Number_of_runs,
                        mc.cores = 1,  # 'mc.cores' > 1 is not supported on Windows
                        function(i) 
                          kmeans(t(Complete_dataaset_220),centers = 4, iter.max = 50,nstart = 5))
# -------------------------------------
# Calculating dissmilarity matrix
K_means_Model <- kmeans(t(Complete_dataaset_220),centers = 4, iter.max = 50,nstart = 10) # nstart: how many random sets should
# Compute all the pairwise dissimilarities (distances) between observations in the data set
Dissimilarity_Matrix <- daisy(t(Complete_dataaset_220),
                              metric = c("euclidean"))# "manhattan", "gower"  
dE2   <- Dissimilarity_Matrix^2
Silhouette_kmeans   <- silhouette(K_means_Model$cluster, dE2)
plot(Silhouette_kmeans)

#----------------------------------------------------------------------------------------------------
# 4) Combine the cluster labels (optional) (next week)
#----------------------------------------------------------------------------------------------------


#####################################################################################################
#                                   Main programme - End
#####################################################################################################