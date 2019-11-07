# NMF analysis
# Dr Reza Rafiee, 29th Oct. 2019
# Analysis: multiple imputation modelling using Amelia and mice packages in R
# The input csv file includes numeric values with some missing data
# The aim is (but not limited) to the following items:
#   1) Read the CSV file including 220 samples with 17 features, 
#   2) running NMF analysis on input dataset and finding the best rank
# All right reserved!
#-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Libraries will be in this section
library(NMF)
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
# Read my input csv file from the working directory
Complete_dataaset_220 <- as.data.frame(read.csv("Complete_Dataset_220_17_CSC3062_RR_2019.csv",row.names = 1))
min(Complete_dataaset_220, na.rm = TRUE) # 0
max(Complete_dataaset_220, na.rm = TRUE) # 0.9994831
class(Complete_dataaset_220)  # [1] "data.frame"
attributes(Complete_dataaset_220)
dim(Complete_dataaset_220) # [1]  17 220
#----------------------------------------------------------------------------------------------------
#                     First analysis: running NMF on complete dataset
#----------------------------------------------------------------------------------------------------
seq1 <- 2    # initial rank
seq2 <- 7    # final rank, this parameter should be carefully seleceted  
noofrun <- 20  # the number of run,  performing 30-50 runs is considered sufficient
# to get a robust estimate of the factorisation rank (Brunet2004; Hutchins2008).
initseed <- 123456
usedmethod <- "ns"
Res_MultiRank <- nmf(Complete_dataaset_220, seq(seq1,seq2), 
                     method=usedmethod, nrun=noofrun, 
                     seed=initseed, .options = "t") # nsNMF is the best algorithm in terms of minimum residual errors
plot(Res_MultiRank)

# which rank is the best? check the NMF rank survey and assess the cophenetic score for different ranks
# The proper factorization rank should be selected where
# the magnitude of the cophenetic correlation coefficient begins to fall.
length(Res_MultiRank$measures$rank) # number of rank used in nmf function
cophen_max <- max(Res_MultiRank$measures$cophenetic)  
for (i in 1:length(Res_MultiRank$measures$rank)) {
  if (Res_MultiRank$measures$cophenetic[i] == cophen_max)
  {
    idx_cophen_max <- i
  }
}  # for loop

# assess the silhouette!
# Res_MultiRank$measures$silhouette.consensus

NMFfitClass <- Res_MultiRank$fit[[idx_cophen_max]]
H_matrix <- NMFfitClass@fit@H
W_matrix <- NMFfitClass@fit@W

# Save your H and W matrix in your working directory
write.csv(H_matrix,"H_matrix_17_220_k4_4_220.csv")
write.csv(W_matrix,"W_matrix_17_220_k4_17_14.csv")

# or you could use the following commands to obtain W and H matrix 
w_matrix <- basis(NMFfitClass)  # basis: W
h_matrix <- coef(NMFfitClass)   # Coef: H
w_matrix %*% h_matrix # ~ mat
dim(w_matrix) # [1] 17  4

# mixture coefficent map for the best rank
coefmap(Res_MultiRank$fit[[idx_cophen_max]])
# mixture coefficent map (heatmap with HC clustering) for the best rank (which is 4)
coefmap(Res_MultiRank$fit[[idx_cophen_max]])

# Consensub map plots for the best rank
consensusmap(Res_MultiRank$fit[[idx_cophen_max]])
filename_NMF_pdf <- paste("NMF_Rank_",idx_cophen_max+1,"Metagenes",nrow(Complete_dataaset_220),"genes_",ncol(Complete_dataaset_220),"samples","_NMF_",cophen_max+1,".pdf",sep = "")
pdf(filename_NMF_pdf)
dev.off()
#par(mfrow=c(2,2))

# If you would like to save your NMF rank survey plot as a PDF file 
#plot(Res_MultiRank)
#filename2_NMF_pdf <- paste("NMF_Consensus_",idx_cophen_max+1,"Metagenes",nrow(Complete_dataaset_220),"genes_",ncol(Complete_dataaset_220),"samples","_NMF_",cophen_max+1,".pdf",sep = "")
#pdf(filename2_NMF_pdf)
#dev.off()


# sample clustering from best fit
plot(silhouette(Res_MultiRank$fit[[idx_cophen_max]]))
par(mfrow=c(2,2))
consensusmap(Res_MultiRank, annCol=NA, labCol=NA, labRow=NA, info=F)

#----------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------- 
# 4. Run PCA algorithm on the factorisation results (i.e., H matrix) obtained from the NMF analysis
# --------------------------------------------------------------------------------------- 
PCA_H_Matrix_Model <- prcomp(t(H_matrix)) # transpose H-matrix (4*220)
summary(PCA_H_Matrix_Model)
PCA_H_Matrix_Model$x

#----------------------------------------------------------------------------------------------------
# 5. Visualise the three main principal component analysis of the H matrix while using 
# colour information in link with the sample subgroup (each sample name includes a subgroup label).
# --------------------------------------------------------------------------------------- 
pca_matrix <- PCA_H_Matrix_Model$x
samples_names <- row.names(PCA_H_Matrix_Model$x)
sample_subgroup_matrix <- matrix(unlist(strsplit(as.character(samples_names), "_\\s*(?=[^_]+$)", perl=TRUE)), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Sample", "Subgroup")))
pca_matrix <- cbind(PCA_H_Matrix_Model$x, "Subgroup"=sample_subgroup_matrix[,2])
pca_matrix[,"Subgroup"][pca_matrix[,"Subgroup"]==1] <- "blue"
pca_matrix[,"Subgroup"][pca_matrix[,"Subgroup"]==2] <- "red"
pca_matrix[,"Subgroup"][pca_matrix[,"Subgroup"]==3] <- "yellow"
pca_matrix[,"Subgroup"][pca_matrix[,"Subgroup"]==4] <- "green"

# 3D plot
pca3d(PCA_H_Matrix_Model$x[,1:3],components = c(1,2,3), col=pca_matrix[,"Subgroup"], group = pca_matrix[,"Subgroup"])

# Pair plot
pairs(PCA_H_Matrix_Model$x, col=pca_matrix[,"Subgroup"])

