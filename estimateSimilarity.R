
#-----Parameter Descriptions------------------------------------#
# - X: If stacked data is provided, the stacked single-cell data for data type 1
# - G: If stacked data is provided, the stacked single-cell data for data type 2
# - y: Vector of outcomes
# - subjectIDs: vector of subjects from 1 to N for the cells
# - outLoc_: If stacked data is NOT provided, the location where the data types are stored for each subject separately
# - zeroAsMissing: Should zeros in the data be treated as missing? T/F
#-------------------------------------------------#

estimateSimilarity <- function(X=NULL,
                               G=NULL, 
                               y=NULL, 
                               subjectIDs=NULL,
                               outLoc_=NULL,
                               zeroAsMissing=FALSE)
  
{
  
  #Sample size
  N <- length(y)
  
  similarityList <- list()
  #Fetch the similarity within each subject between X and G#Initialize
  for(k in 1:N){
    #Initialize
    if(!is.null(outLoc_)){
      X_k <- readRDS(paste0(outLoc_, "/X_", k, ".rda"))
      G_k <- readRDS(paste0(outLoc_, "/G_", k, ".rda"))
    }else {
      X_k <- X[subjectIDs==k,]
      G_k <- G[subjectIDs==k,]
    }
    if(zeroAsMissing){
      X_k[X_k==0] <- NA
      G_k[G_k==0] <- NA
    }
    #Calculate the similarity
    similarityList[[k]] <- cor(x=X_k, y=G_k, method="pearson", use="pairwise.complete.obs")
    #Convert any NAs to 0s
    similarityList[[k]][is.na(similarityList[[k]])] <- 0
  }
  #Find how many features in X and G
  p <- ncol(X_k)
  q <- ncol(G_k)
  
  #Put the similarity matrix list into a tensor
  similarityArray <- array(data=unlist(similarityList), dim=c(p,q,N))
  Z_tensor <- rTensor::as.tensor(similarityArray)
  
  return(Z_tensor)
  
}