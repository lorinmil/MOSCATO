

# Z_tensor: Similarity tensor to be used as predictor tensor
# y: Vector of outcomes (one per subject)
# distY: Distribution of y to be used in GLM
# numIterForStability: Number of subsamples to be used for checking stability
# sizeOfSubsampleForStability: OPTIONAL. Number of features to use for each subsample. If empty, defaults to half of the subjects.
# phi: Instability threshold
# minSelect: Minimum number of features to be selected. Vector of size 2 where the first number is for X and the other is for G
# epsilon: Difference for establishing convergence
# minIter: Minimum number of iterations to do before considering convergence
# maxIter: Maximum number of iterations before deciding it didn't converge
# numAlpha: OPTIONAL: Number of alpha's to consider between 0 and 1. But, if alphaGrid is specified, that will be used instead
# alphaGrid: Vector of alpha values to consider for tuning

tuneMOSCATO <- function(#Inputs for data
                        Z_tensor, #Similarity tensor
                        y, #Vector of outcomes
                        distY="gaussian", #Distribution of y
                        #Inputs for tuning using StARS
                        numIterForStability=20, #Number of subsamples to establish stability
                        sizeOfSubsampleForStability=NULL, #The number of subjects for each subsample
                        phi=0.05, #Instability threshold
                        minSelect=NULL, #Minimum number of variables to select
                        epsilon=0.1, #Convergence criteria
                        minIter=5, #Minimum number of iterations when fitting tensor regression
                        maxIter=50, #Maximum number of iterations when fitting tensor regression)
                        numAlpha = 5, #How many alpha values to check for between 0 and 1 in the tuning?
                        alphaGrid = NULL #Specifically which alpha to run?
){
  
  
  #Number of observations
  N <- length(y)
  
  #----------------Pre-specify the subsamples-----------------#
  subsample <- list()
  if(is.null(sizeOfSubsampleForStability)){
    sizeOfSubsampleForStability <- floor(0.5*N)
  }
  if(choose(N, sizeOfSubsampleForStability) <= numIterForStability){
    allCombos <- combn(N, sizeOfSubsampleForStability)
    subsample <- list()
    for(iter in 1:ncol(allCombos)){
      subsample[[iter]] <- allCombos[,iter]
    }
  }else {
    for(iter in 1:numIterForStability){
      
      #Subset the data to a random sample
      subsample[[iter]] <- base::sample(1:N, size=sizeOfSubsampleForStability)
    }
  }
  
  
  
  
  
  
  #----------------Tune values for tensor regression------------------#
  
  #Initialize the values
  D <- length(dim(Z_tensor))-1
  dims <- dim(Z_tensor)[1:D]
  maxSelect <- dims
  if(is.null(minSelect)){
    minSelect <- rep(1, D)
  }
  alphaVals <- rep(0.5, D)
  alphaList <- list()
  jList <- list()
  DhatList <- list()
  
  
  #Loop through the maximum values until optimal instability
  for(d in 1:D){
    Dhat <- c()
    alpha <- c()
    jVec <- c()
    for(j in minSelect[d]:dims[d]){
      
      cat("Processing dimension ", d, "with max vars set to ", j, "...\n")
      
      #Update
      maxSelect[d] <- j
      
      #---Tune alpha for this number of maximum variables---#
      tuned <- function_tuneAlpha(X_=Z_tensor,
                                  y_=y,
                                  distY=distY,
                                  subsample_=subsample,
                                  phi_=phi,
                                  maxSelect_=maxSelect,
                                  alphaVals_=alphaVals,
                                  dim_=d,
                                  epsilon=epsilon,
                                  minIter=minIter,
                                  maxIter=maxIter,
                                  numAlpha = numAlpha,
                                  alphaGrid = alphaGrid)
      
      if(tuned$alphaConverged){
        cat("Alpha converged! Dhat=", tuned$Dhat, "alpha=", tuned$selectAlpha, "max=", j, "/n")
        Dhat <- c(Dhat, tuned$Dhat)
        alpha <- c(alpha, tuned$selectAlpha)
        jVec <- c(jVec, j)
        
        #Update the list
        alphaList[[d]] <- alpha
        jList[[d]] <- jVec
        DhatList[[d]] <- Dhat
        
        #If this is greater than phi, then 
        if(Dhat[length(Dhat)]>phi){
          if(length(Dhat)>1){
            cat("Maximum selected!")
            #Update the alpha with the tuned alpha value
            alphaVals[d] <- alpha[length(Dhat)-1]
            maxSelect[d] <- jVec[length(Dhat)-1]
            break
          }else {
            cat("Maximum selected!")
            #Update the alpha with the tuned alpha value
            alphaVals[d] <- alpha[length(Dhat)]
            maxSelect[d] <- jVec[length(Dhat)]
            break
          }
        }
      } #End if alpha converged
      
    } #Loop to the next maximum number of features
  } #Loop to the next dimension
  
  
  
  
  #Compile objects to return from the function
  tuningList <- list(alpha = alphaList,
                     Dhat = DhatList,
                     max = jList,
                     phi = phi,
                     subsample = subsample,
                     numAlpha = numAlpha,
                     alphaGrid = alphaGrid,
                     tunedAlpha = alphaVals,
                     tunedMax = maxSelect)
  
  return(tuningList)
}