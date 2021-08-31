

#---------------Purpose: Tune the alpha from the tensor elastic net model-----------#

function_tuneAlpha <- function(X_,         #Predictor data
                               y_,         #Outcome
                               distY, #Distribution for the outcome
                               subsample_, #List of subsamples to make
                               phi_,      #phi to use for instability
                               maxSelect_,  #Vector of maximum values to use for the dimensions
                               alphaVals_,  #Vector of alphas to use (this is only used for the dimensions outside of dim_)
                               dim_,         #Dimension to tune alpha for
                               #Parameters for the tensor regression
                               epsilon,
                               minIter,
                               maxIter,
                               numAlpha = 10,
                               alphaGrid = NULL
){
  
  #Initialize
  Dhat <- c() #This is the set of instability that has been processed for corresponding processed alphas
  alpha_use <- c() #This is the set of alpha values that have been checked
  if(is.null(alphaGrid)){
    alpha_vec <- seq(from=0.01, to=1, by=1/numAlpha) #This is the cumulative set of alphas to check
  } else{
    alpha_vec <- alphaGrid
  }
  alpha_vec_tmp <- c()  #This is the set of alphas that have been processed
  glmnetAlpha <- alphaVals_ #Values to use for alpha across all dimensions (only the ones corresponding outside dim_ will be used)
  instabilityDF <- NULL #Just inialize in case nothing converges
  
  #How many subsamples?
  numSubsamples <- length(subsample_)
  
  #How many variables on each dimension?
  D <- length(dim(X_))-1
  dims <- dim(X_)[1:D]
  
  #Intialize
  alpha_vec_toAdd <- c()
  
  #Loop through the alpha values
  # alphaResults <- lapply(alpha_vec, FUN=function(alpha){
  for(alpha in alpha_vec){
    
    cat(paste("Processing alpha: ", alpha, "..."))
    
    #Update which alpha to process
    glmnetAlpha[dim_] <- alpha
    
    #Loop through all of the subsamples
    subsampleResults <- mclapply(X=1:numSubsamples, FUN=function(sub){
      #Subset the data to a random sample
      subsample_iter <- subsample_[[sub]]
      Z_tensor_iter <- X_[,,subsample_iter]
      y_iter <- y_[subsample_iter]
      
      #Perform the tensor regression on the subsample
      successRun <- tryCatch({
        tensorReg <- tensorElasticNet(Z_tensor=Z_tensor_iter,
                                      y=y_iter,
                                      distY=distY,
                                      epsilon=epsilon,
                                      minIter=minIter,
                                      maxIter=maxIter,
                                      glmnetAlpha=glmnetAlpha,
                                      maxSelect=maxSelect_)
        modelConverged <- tensorReg$fittedResults$converged
        erroredOut <- FALSE
      }, error=function(cond){ #If there is an error in performing glmnet
        cat("Error!")
        return("Failure")
      })
      
      if(successRun=="Failure"){
        erroredOut <- TRUE
        modelConverged <- FALSE
      }
      
      #Only process if it converged
      if(!is.null(modelConverged) && !is.na(modelConverged) && modelConverged){
        
        #Mark which were selected with 0s and 1s
        u <- rep(0, length(tensorReg$finalBeta[[dim_]]))
        u[tensorReg$finalBeta[[dim_]]!=0] <- 1
        
      }else { #End if the model converged
        modelConverged <- FALSE
        #Leave this here as a dummy vector (only used when errors out)
        u <- rep(0, dims[dim_])
      }
      
      return(list(u=u, modelConverged=modelConverged, erroredOut=erroredOut))
    }) #Loop to the next subsample
    
    #How many converged?
    converged <- c()
    for(sub in 1:numSubsamples){
      converged <- c(converged, subsampleResults[[sub]][[2]])
    }
    numConverged <- sum(converged)
    
    #Only calculate instability if enough runs were made
    if(numConverged>=floor(0.75*numSubsamples)){
      
      #Compile subsamples
      first <- TRUE
      for(sub in 1:numSubsamples){
        if(subsampleResults[[sub]][[2]] || subsampleResults[[sub]][[3]]){
          if(first){
            first <- FALSE
            selected <- subsampleResults[[sub]][[1]]
          }else {
            selected <- rbind(selected, subsampleResults[[sub]][[1]])
          }
        }
      } #Loop to next subsample
      
      #Calculate the average instability
      ColMeans <- colMeans(selected)
      xis <- 2 * ColMeans * (1 - ColMeans)
      Dhat <- c(Dhat, mean(xis))
      alpha_use <- c(alpha_use, alpha)
    }
    
  } #End loop for alpha list
  
  
  
  #If no alphas converged, then mark as such
  if(is.null(Dhat)){
    alphaConverged <- FALSE
    selectAlpha <- NULL
    selectDhat <- NULL
  }else { #Otherwise, find the best alpha
    alphaConverged <- TRUE
    if(sum(Dhat<phi_)>0){
      selectAlpha <- alpha_use[which.min(Dhat)]
      selectDhat <- Dhat[which.min(Dhat)]
    }else {
      #Otherwise take the smallest instability if nothing is less than phi
      selectAlpha <- alpha_use[which.min(Dhat)]
      selectDhat <- Dhat[which.min(Dhat)]
    }
  }
  
  
  #Return the values
  returnList <- list(alphaConverged=alphaConverged,
                     selectAlpha=selectAlpha,
                     Dhat=selectDhat)
  return(returnList)
  
}












#---------------Purpose: Perform elastic net tensor regression-----------#

tensorElasticNet <- function(Z_tensor,     #Tensor for the independent variables
                             y,           #Outcome vector
                             distY="gaussian", #Distrubition for the outcome
                             epsilon=0.1,     #Convergence criteria
                             minIter=5,     #Minimum number of iterations
                             maxIter=50,     #Maximum number of iterations
                             glmnetAlpha, #The alpha parameter for glmnet (the weight put on lasso)
                             glmnetLambda, #The lambda parameters for glmnet (the sparsity parameters)
                             maxSelect=NULL #Maximum number of features to select (if this is specified, lambda is ignored)
){
  
  #Initialize the values
  D <- length(dim(Z_tensor))-1
  dims <- dim(Z_tensor)[1:D]
  p <- dims[1]
  q <- dims[2]
  betaList <- list()
  alphaList <- list()
  lambda <- rep(0, D)
  
  #Initialize to zero
  numIter <- 0
  
  #Set the random intialization
  beta_t <- list()
  for(d in 1:D){
    beta_t[[d]] <- rnorm(dims[d], 0, 1)
  }
  
  #Initialize
  alpha_t <- as.numeric(glm(y ~ 1, family=distY)$coefficients)
  
  #Loop until convergence
  repeat{
    #Initialize for the loop
    beta_t2 <- beta_t
    #Iteratively update the dimensions
    for(d in 1:D){
      #Multiply out the other dimensions
      X_prime <- FUNC_gradient_wk_fWb_X(W_=unlist(beta_t2),
                                        X_=Z_tensor,
                                        k_=d)
      colnames(X_prime) <- paste0("x", 1:ncol(X_prime))
      #Execute the final model
      if(is.null(maxSelect)){
        model_t <- glmnet::glmnet(y=y, 
                                  offset=rep(alpha_t, length(y)),
                                  x=X_prime,
                                  family=distY,
                                  alpha=glmnetAlpha[d],
                                  standardize=FALSE,
                                  intercept=FALSE,
                                  lambda = glmnetLambda[d])
        beta_t2[[d]] <- as.numeric(model_t$beta)
        beta_t2[[d]][is.na(beta_t2[[d]])] <- beta_t[[d]][is.na(beta_t2[[d]])]
      }else {
        model_t <- glmnet::glmnet(y=y, 
                                  offset=rep(alpha_t, length(y)),
                                  x=X_prime,
                                  family=distY,
                                  alpha=glmnetAlpha[d],
                                  standardize=FALSE,
                                  intercept=FALSE,
                                  dfmax=maxSelect[d])
        beta_t2[[d]] <- as.numeric(model_t$beta[,max(which(model_t$df<=maxSelect[d]))])
        beta_t2[[d]][is.na(beta_t2[[d]])] <- beta_t[[d]][is.na(beta_t2[[d]])]
        lambda[d] <- model_t$lambda[max(which(model_t$df<=maxSelect[d]))]
      }
    } #Loop to the next dimension
    
    ###Update alpha
    #Multiply out all coefficients
    X_prime <- FUNC_f_Wb_X(W_=unlist(beta_t2), X_=Z_tensor)
    alpha_t <- glm(y ~ 1, 
                   offset=X_prime,
                   family="gaussian")$coefficients[1]
    
    ###Check for convergence
    #New likelihood
    logLik2 <- log(length(y)/sqrt(2*pi))-0.5*sum((y-(alpha_t+X_prime))^2)
    
    #Add the values to the list to check for convergence afterwards
    if(numIter==0){
      logLikeList <- logLik2
      alphaList <- alpha_t
      betaList <- matrix(unlist(beta_t), nrow=1, ncol=length(unlist(beta_t)))
    }else {
      logLikeList <- c(logLikeList, logLik2)
      alphaList <- c(alphaList, alpha_t)
      betaList <- rbind(betaList, unlist(beta_t))
    }
    
    #Check if the likelihood has converged
    if(numIter>minIter){
      if((logLik2>logLik1) && (abs(logLik2 - logLik1) < epsilon)){
        finalBeta <- beta_t2
        finalAlpha <- alpha_t
        fittedResults <- data.frame(converged=TRUE, likelihood=logLik2, numIterations=numIter, alpha=alpha_t)
        cat("Converged!\n")
        break
      }
    }
    
    #Break out if the maximum number of loops is hit
    if(numIter >= maxIter){
      #Take the mean of last iterations if it did not converge
      beta_t_new <- colMeans(betaList[(nrow(betaList)-10):nrow(betaList),])
      beta_t2 <- list()
      startIdx <- 1
      for(d in 1:D){
        beta_t2[[d]] <- beta_t_new[startIdx:(startIdx + dims[d] - 1)]
        startIdx <- startIdx + dims[d]
      }
      finalBeta <- beta_t2
      finalAlpha <- alpha_t
      fittedResults <- data.frame(converged=FALSE, likelihood=logLik2, numIterations=numIter, alpha=alpha_t)
      cat("Did not converge!\n")
      break
    }
    
    #Increment values
    numIter <- numIter + 1
    beta_t <- beta_t2
    logLik1 <- logLik2
    
  } #Loop to the next interation for convergence
  
  #Mark which variables were selected
  returnX <- rep(0, length(finalBeta[[1]]))
  returnG <- rep(0, length(finalBeta[[2]]))
  returnX[which(finalBeta[[1]]!=0)] <- 1
  returnG[which(finalBeta[[2]]!=0)] <- 1
  
  
  #Return the values
  return(list(returnX=returnX,
              returnG=returnG,
              fittedResults=fittedResults,
              finalBeta=finalBeta,
              betaList=betaList,
              logLikeList=logLikeList,
              lambda = lambda))
  
  
} #End the tensorElasticNet function













#Function for calculating the gradient_wk of f_Wb(Xi). So multiply out all loadings from X except for the dimension k_ for each subject
#INPUTS:  W_  - Vector of loadings
#         X_  - Input tensor for predictors
#         k_  - The dimension for the gradient  
#RETURN:  Matrix of dimensions N by number of variables in dimension k_
FUNC_gradient_wk_fWb_X <- function(W_, X_, k_){
  
  #Fetch the dimensions of X
  dims <- dim(X_)
  
  #Fetch the sample size
  N <- dim(X_)[length(dim(X_))]
  
  #Fetch the number of dimensions (excluding the subject dimension) in the input tensor
  K <- length(dim(X_)) - 1
  
  #Initialize gradient matrix
  gradient_wk_fWb_X <- matrix(rep(0,N*dims[k_]), nrow=N, ncol=dims[k_])
  
  #Loop through each subject
  # gradient_wk_fWb_X <- mclapply(X=1:N, #mc.cores = detectCores()/2-1, mc.preschedule=FALSE, 
  #                               FUN=function(i){
  for(i in 1:N){
    #Fetch the tensor for subject i (TODO: make this more robust to higher dimensions)
    if(K==2){
      gradient_wk_fWb_X_i <- X_[,,i]
    } else if(K==3){
      gradient_wk_fWb_X_i <- X_[,,,i]
    } else if(K==4){
      gradient_wk_fWb_X_i <- X_[,,,,i]
    }
    #Calculate the partial gradiant after removing the kth coeff vector
    for(k_2 in 1:K){
      if(k_2==1){
        idx_k2 <- 1:dims[k_2]
      } else {
        idx_k2 <- (sum(dims[1:(k_2-1)]) + 1):(sum(dims[1:(k_2-1)]) + dims[k_2])
      }
      #Don't multiply the input dimension k_
      if(k_2 != k_){
        w_k2 <- as.matrix(W_[idx_k2])
        #Check to see if the w vector should take the transpose or not
        if(dim(gradient_wk_fWb_X_i)[k_2]==nrow(w_k2)){
          gradient_wk_fWb_X_i <- rTensor::ttm(tnsr=gradient_wk_fWb_X_i,mat=t(w_k2), m=k_2)
        }else {
          gradient_wk_fWb_X_i <- rTensor::ttm(tnsr=gradient_wk_fWb_X_i,mat=w_k2, m=k_2)
        }
      }
    } #Loop to the next k2 dimension
    #Turn the gradient into a vector instead of a tensor and add it to the vector of values
    # return(vec(gradient_wk_fWb_X_i))
    
    gradient_wk_fWb_X[i,] <- vec(gradient_wk_fWb_X_i)
  }
  # )
  
  #Put into a matrix 
  # gradient_wk_fWb_X <- matrix(unlist(gradient_wk_fWb_X), nrow=N, ncol=dims[k_], byrow=TRUE)
  
  return(gradient_wk_fWb_X)
}













#Function for calculating f_Wb(X). Multiplies out the loadings 
#INPUTS:  X_  - Input tensor
#         W_  - vector of loadings
#RETURNS: Vector of size N of the subject f_Wb(Xi) values
FUNC_f_Wb_X <- function(X_, W_){
  
  #Fetch the dimensions of X
  dims <- dim(X_)
  
  #Fetch the sample size
  N <- dim(X_)[length(dim(X_))]
  
  #Fetch the number of dimensions (excluding the subject dimension) in the input tensor
  K <- length(dim(X_)) - 1
  
  #Multiply X by all of the loadings to get f(W,b)(Xi) for each subject i
  f_wbX <- rep(0,N) #initialize
  
  # f_wbX <- mclapply(X=1:N, mc.cores = detectCores()/2-1, mc.preschedule=FALSE, FUN=function(i){
  for(i in 1:N){
    #Fetch the tensor for subject i (TODO: make this more robust to higher dimensions)
    if(K==2){
      f_wbX_i <- X_[,,i]
    } else if(K==3){
      f_wbX_i <- X_[,,,i]
    } else if(K==4){
      f_wbX_i <- X_[,,,,i]
    }
    #Loop through the dimensions and multiply by the loading
    for(k_2 in 1:K){
      if(k_2==1){
        idx_k2 <- 1:dims[k_2]
      } else {
        idx_k2 <- (sum(dims[1:(k_2-1)]) + 1):(sum(dims[1:(k_2-1)]) + dims[k_2])
      }
      #Use the current loadings if they've already been processed
      w_tminus1_k2 <- as.matrix(W_[idx_k2])
      #Check to see if the w vector should take the transpose or not
      if(dim(f_wbX_i)[k_2]==nrow(w_tminus1_k2)){
        f_wbX_i <- rTensor::ttm(tnsr=f_wbX_i,mat=t(w_tminus1_k2), m=k_2)
      }else {
        f_wbX_i <- rTensor::ttm(tnsr=f_wbX_i,mat=w_tminus1_k2, m=k_2)
      }
    } #Loop to the next k2 dimension
    
    #Turn the gradient into a vector instead of a tensor and add it to the vector of values
    # return(vec(f_wbX_i))
    f_wbX[i] <- vec(f_wbX_i)
  }
  # )
  
  # f_wbX <- unlist(f_wbX)
  
  return(f_wbX)
  
}

