

#-----------------------------------------------------#
#--PURPOSE: Simulate Single Cell Sparse Networks------#
#-----------------------------------------------------#


####INPUTS FOR SPLATTER:###
####Gene mean###
# pi0       - Outlier Probability
# mu0       - outlier location
# sigma0    - Outlier scale
####Cell Mean###
# mu_L      - Library size location
# sigma_L   - Library size scale
####Trended Cell Mean
# phi       - common dispersion
# df0       - BCV DoF
####Final Count
# x0        - dropout modpoint
# k_shape   - Dropout shape
####Sample-level inputs
# m_mean    - Mean number of cells per subject
############################



####INPUTS FOR FANS:####
# alpha     - Helps specify strength of network
# Rsq       - Helps specify strength of network
# N         - Number of subjects to simulate
# rind      - Helps specify strength of network
# nx        - Number of network features within X
# nxPrime   - Number of features within X that are related to G but not the outcome
# ns        - Number of features within X that are related to the outcome but not G
# nxNoise   - Number of features within X that are not related to G or the outcome
# ng        - Number of network features within G
# ngPrime   - Number of features within G that are related to X but not the outcome
# nh        - Number of features within G that are related to the outcome but not X
# ngNoise   - Number of features within G that are not related to the outcome or X
# noiseSD   - The standard deviation for the noise features
# noiseMean - The mean for the noise features
########################





#---------------------------------------------------------------------------------------#
#--------------Define the functionfor simulating the single cell network----------------#
#---------------------------------------------------------------------------------------#

simSingleCellNetwork <- function(#Inputs for Splatter
                                 pi0, mu0, sigma0, mu_L, sigma_L, phi, df0, x0, k_shape, m_mean,
                                 #Inputs for Network
                                 alpha, Rsq, N, rind, nx, nxPrime, ns, nxNoise, ng, ngPrime, nh, ngNoise, noiseSD, noiseMean,
                                 outLoc_=NULL
){
  
  #------Set up covariance structure for latentvariables-------#
  tauInd <- (1 - rind^2)/rind^2
  
  #Set  up the covariance matrix for (Y, X, G, S, G, X', G')
  sigma <- calSigma(setPar(alpha, Rsq)) #Set up the standard devs for the latent vars
  
  
  
  
  #----------------Initialize objects for simulation-----------#
  
  #Find the total number of features in X and G (p and q, respectively)
  p <- nx + nxPrime + ns + nxNoise
  q <- ng + ngPrime + nh + ngNoise
  
  
  
  #----------------Simulate values for each subject------------#
  
  #Initilalize
  outcome <- c()
  subjectIDs <- c()
  
  #Loop through and generate values for each subject
  for(i in 1:N){
      
    
    #-----------Get gene mean profiles------------#
    #Step 0: Set up latent variables
    latent <- (genLatent(sigma)/sqrt(diag(sigma)))^2
    
    ##SIMULATE Y##
    y <- rnorm(1, latent[1], sqrt(tauInd * sigma["Yprime", "Yprime"]))
    
    ##SIMULATE X GENE MEANS##
    #Step 1: Original Mean
    lamda_prime_X <- c(rep(latent[2], nx),
                       rep(latent[6], nxPrime),
                       rep(latent[4], ns),
                       rnorm(n=nxNoise, mean=noiseMean, sd=noiseSD)^2)
    
    #Step 2: Gene Mean
    #Outlier indicator
    one0_X <- rbinom(p, 1, pi0)
    #Outlier factor
    phi_prime_X <- rlnorm(p, mu0, sigma0)
    
    lamda_X <- c()
    for(k in 1:p){
      #Gene mean
      lamda_X <- c(lamda_X, one0_X[k]*(phi_prime_X[k]*median(lamda_prime_X)) + (1-one0_X[k])*lamda_prime_X[k])
    }
    
    ##SIMULATE G GENE MEANS##
    #Step 1: Original Mean
    lamda_prime_G <- c(rep(latent[3], ng),
                       rep(latent[7], ngPrime),
                       rep(latent[5], nh),
                       rnorm(n=ngNoise, mean=noiseMean, sd=noiseSD)^2)
    
    #Step 2: Gene Mean
    #Outlier indicator
    one0_G <- rbinom(q, 1, pi0)
    #Outlier factor
    phi_prime_G <- rlnorm(q, mu0, sigma0)
    
    lamda_G <- c()
    for(k in 1:q){
      #Gene mean
      lamda_G <- c(lamda_G, one0_G[k]*(phi_prime_G[k]*median(lamda_prime_G)) + (1-one0_G[k])*lamda_prime_G[k])
    }
    
    
    #Combine the X and G lamdas
    lamda <- c(lamda_X, lamda_G)
    
    
    #-----------Simulate the actual values------------#
    #Generate how many cells they should have
    m <- rpois(1, m_mean)
    
    #Initialize objects
    subj <- rep(i, m)
    
    XG_i <- matrix(nrow=m, ncol=p+q+1)
    colnames(XG_i) <- c("subject", c(paste0("x",1:p), paste0("g",1:q)))
    XG_prime_i <- matrix(nrow = m, ncol = p+q+1)
    colnames(XG_prime_i) <- c("subject", c(paste0("x",1:p), paste0("g",1:q)))
    XG_i[,1] <- subj
    XG_prime_i[,1] <- subj
    
    BCV <- rchisq(n=1, df=df0)
    
    #Loop through to find values in every cell
    for(j in 1:m){
      
      #Step 3: Cell mean
      #Expected library size
      L_j <- rlnorm(1, mu_L, sigma_L)
      
      lamda_prime_ij <- L_j*lamda/sum(lamda)
      
      B_ij <- (phi + 1/sqrt(lamda_prime_ij))*sqrt(df0/BCV)
      
      B_lamda <- cbind(B_ij, lamda_prime_ij)
      B_lamda <- unique(B_lamda)
      
      lamda_ij <- c()
      for(idx in 1:nrow(B_lamda)){
        lamda_ij <- c(lamda_ij, rep(rgamma(n=1, 
                                           shape=(1/(B_lamda[idx,1]^2)), 
                                           scale=B_lamda[idx,2]*(B_lamda[idx,1]^2)), 
                                    sum(B_ij==B_lamda[idx,1] & lamda_prime_ij==B_lamda[idx,2])))
      }
      
      #Simulate the raw values of X in cell j
      for(k in 1:(p+q)){
        
        #Step 5: True Count
        XG_prime_i[j,k+1] <- rpois(1, lamda_ij[k])
        
        #Step 6: Final Count
        #Dropout probability
        pi_D_ij <- 1/(1+exp(-k_shape*(log(lamda_ij[k])-x0)))
        #Dropout Indicator
        # one_D_ij <- rbinom(1, 1, 1-pi_D_ij)
        one_D_ij <- rbinom(1, 1, pi_D_ij)
        #Count
        XG_i[j,k+1] <- one_D_ij*XG_prime_i[j,k+1] 
        
      } #Loop to the next feature
      
    } #Loop to the next cell in subject i
    
    #If the output directory is specified, save the simulation in its own file
    if(!is.null(outLoc_)){
      saveRDS(XG_i[,2:(p+1)], paste0(outLoc_, "/X_", i, ".rda"))
      saveRDS(XG_i[,c((p+2):ncol(XG_i))], paste0(outLoc_, "/G_", i, ".rda"))
    }else {
      subjectIDs <- c(subjectIDs, XG_i[,1])
      if(i==1){
        X <- XG_i[,2:(p+1)]
        X_prime <- XG_prime_i[,2:(p+1)]
        
        G <- XG_i[,(p+2):ncol(XG_i)]
        G_prime <- XG_prime_i[,(p+2):ncol(XG_i)]
      }else {
        X <- rbind(X, XG_i[,2:(p+1)])
        X_prime <- rbind(X_prime, XG_prime_i[,2:(p+1)])
        
        G <- rbind(G, XG_i[,(p+2):ncol(XG_i)])
        G_prime <- rbind(G_prime, XG_prime_i[,(p+2):ncol(XG_i)])
      }
    }
    
    outcome <- c(outcome, y)
    
    #Log message
    cat("Finished Subject: ", i, "!\n")
    
  }
  
  #Save it
  y <- outcome
  
  
  
  
  #----------------------Return the simulated values-----------------#
  xDescription <- c(rep("Network", nx), rep("xPrime", nxPrime), rep("S", ns), rep("Noise", nxNoise))
  gDescription <- c(rep("Network", ng), rep("gPrime", ngPrime), rep("H", nh), rep("Noise", ngNoise))
  #Save the other info
  if(!is.null(outLoc_)){
    saveRDS(y, paste0(outLoc_, "/Ysim.rda"))
    saveRDS(list(xDescription=xDescription, gDescription=gDescription), paste0(outLoc_, "/description.rda"))
  } else{
    return(list(G=G, 
                X=X, 
                subjectIDs=subjectIDs,
                y=y,
                xDescription=xDescription, 
                gDescription=gDescription))
    
  }
  
} #End the simSingleCellNetwork function

