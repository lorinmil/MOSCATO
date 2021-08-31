



calSigma <- function(sysPar){
  
  cc <- unlist(sysPar)
  for( i in 1:length(cc) ){ 
    assign(names(cc)[i],cc[i]) 
  }
  
  #Covariance matrix for (Y, X, G, S, H)
  Sigmaprime <- matrix(c(syxgsh + sxg * alpha + sg + alpha^2 * ss + alpha^2 * sh, # 1,1
                         sxg * sqrt(alpha) + sg * sqrt(alpha), # 1,2
                         sg, # 1,3
                         alpha * ss, #1,4
                         alpha * sh, #1,5
                         sxg * sqrt(alpha) + sg * sqrt(alpha), #2,1
                         sxg + sg * alpha, #2,2
                         sg * sqrt(alpha), #2,3
                         0, #2,4
                         0, #2,5
                         sg, #3,1
                         sg * sqrt(alpha), #3,2
                         sg, #3,3
                         0, #3,4
                         0, #3,5
                         alpha * ss, #4,1
                         0, #4,2
                         0, #4,3
                         ss, #4,4
                         0, #4,5
                         alpha * sh, #5,1
                         0, #5,2
                         0, #5,3
                         0, #5,4
                         sh ), 5,5)
  SigInd <- Sigmaprime[2:3, 2:3] #Elements from covariance matrix for S and H
  Sig2Prime <- matrix(0, nrow=7, ncol=7) #Initialize cov matrix for (Y,X,G,S,H,X',G')
  Sig2Prime[1:5, 1:5] <- Sigmaprime
  Sig2Prime[6:7, 6:7] <- SigInd
  rownames(Sig2Prime) <- c("Yprime", "X", "G", "S", "H", "Xprime", "Gprime")
  colnames(Sig2Prime) <- c("Yprime", "X", "G", "S", "H", "Xprime", "Gprime")
  return(Sig2Prime)
}






setPar <- function(alpha, Rsq)
{
  sg <- 1  # var(G) set to 1
  
  sh <- sg    # same as G. H is assoc with Y, ow independent
  bxg <- sqrt(alpha)
  byx.g <- sqrt(alpha)
  byg.x <- 1 - alpha
  bys <- alpha
  byh <- alpha
  tau <- (1 - Rsq)/Rsq
  sxg <- alpha * tau   # error variance of X given G
  ss <- sxg + sg * alpha    # same as X. S is assoc with Y, ow independent
  
  syxgsh <- tau + (alpha * tau)^2 + 2 * (alpha ^ 2) * tau #alpha^2*tau^2 + 3*tau  # error variance of Y given X, S, G, H
  list(alpha=alpha, Rsq=Rsq, sg=sg, ss=ss, sh=sh, 
       bxg=bxg, 	byx.g=byx.g, byg.x=byg.x, bys=bys, byh=byh, 
       tau=tau, sxg=sxg, syxgsh=syxgsh)
}






genLatent <- function(sigma2prime){
  
  latent <- MASS::mvrnorm(1, rep(0, 7), sigma2prime)
  
  return(latent)
}