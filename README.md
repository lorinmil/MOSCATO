# MOSCATO

## Introduction

This code executes the Multi-Omic Single-Cell Analysis using TensOr regression (MOSCATO) method which uses multi-omic, single-cell data with two data types, X and G, and a univariate outcome which follows an exponential family distribution. In summary, it will estimate the similarity and perform tensor regression using the similarity as the predictor tensor. The features will be selected by applying an elastic net constraint on the coefficient vectors, and any non-zero coefficients will correspond to the features selected. In addition to applying MOSCATO, code to simulate multi-omic, single-cell data is also provided.

## STEP 1: Setting up the data

To use MOSCATO, two data types and a univariate outcome are assumed for each subject. The data types may be supplied either by each subject containing it's own datasets or in the long format where the first column specifies the subject. Note that the data types must have their rows matched by the same cells (i.e., X and G must have the same number of rows where the rows correspond to cells).

If wanting to use simulated data, use the simSingleCellNetwork() function from the simSingleCellNetwork.R file. This will use Splatter simulations proposed by Zappia et al in 2017 with extensions using concepts from the proposal of DNSMI by Zhang et al in 2020. To perform the function, a series of inputs will be specified which correspond to Splatter inputs for single-cell data:

- pi0: Probability of the expression being an outlier?
- mu0: If the expression is an outlier, what is the location?
- sigma0: If the expression is an outlier, what is the scale?
- mu_L: Location for library size
- sigma_L: Scale for library size
- phi: Common dispersion
- df0: BCV degrees of freedom
- x0: Midpoint for determining if dropout (i.e., expression truncated to 0)
- k_shape: Shape for determining if dropout
- m_mean: Mean number of cells for a given subject

There are also parameters for specifying the network strength and size of the simulations, similarly as used in DNSMI:

- alpha: Helps specify signal strength. A number between 0 and 1
- Rsq: Helps specify signal strength. A number between 0 and 1
- rind: Helps specify signal strength. A number between 0 and 1
- N: Number of subjects
- nx: Number of network features in X (i.e., related to features in G and to the outcome)
- nxPrime: Number of features in X related to features within G but not the outcome
- ns: Number of features in X related to the outcome but not any features in G
- nxNoise: Number of features in X unrelated to G and the outcome
- ng: Number of network features in G (i.e., related to features in X and to the outcome)
- ngPrime: Number of features in G related to features within X but not the outcome
- nh: Number of features in G related to the outcome but not any features in X
- ngNoise: Number of features in G unrelated to X and the outcome
- noiseSD: Standard deviation to use for noise features
- noiseMean: Mean to use for the noise features

One other parameter outLoc_ may be specified as a folder to save the data if the subject-level data is wanting to be saved as separate files. Note that this will save files such as X_1.rda and G_1.rda for subject 1, X_2.rda and G_2.rda for subject 2, and so on.

An example is as follows to generate multi-omic single-cell data with a network containing subset of features. 

OPTION 1: Save all subject-level data in separate files

```
  source("supportingSimulationFunctions.R")   
  source("simSingleCellNetwork.R")
  
  outputLocation <- "OPTIONAL: Put output location here for subject files"
  
  testSim <- simSingleCellNetwork(
    #Splatter inputs
    pi0=0.002, #Outlier Probability
    mu0=5, #outlier location
    sigma0=0.4, #Outlier scale
    mu_L=12, #Library size location
    sigma_L=0.2, #Library size scale
    phi=0.1, #common dispersion
    df0=7, #BCV DoF
    x0=1, #dropout modpoint
    k_shape=0.5, #Dropout shape
    #Sample-level inputs
    m_mean=500, #Mean number of cells per subject
    #DNSMI Inputs
    alpha=0.35, 
    Rsq=0.85,
    N=50,
    rind=0.85,
    nx=15,
    nxPrime=20,
    ns=20,
    nxNoise=1500,
    ng=10,
    ngPrime=15,
    nh=15,
    ngNoise=1400,
    noiseSD=1,
    noiseMean=0,
    outLoc_=outputLocation
  )
```

If outLoc_ was specified, this will create files for each of the subjects for their data types, along with a Ysim.rda file for a vector of their continuous outcomes.


OPTION 2: Stack the single-cell data:

```
  source("supportingSimulationFunctions.R")  
  source("simSingleCellNetwork.R")
  
  testSim <- simSingleCellNetwork(
    #Splatter inputs
    pi0=0.002, #Outlier Probability
    mu0=5, #outlier location
    sigma0=0.4, #Outlier scale
    mu_L=12, #Library size location
    sigma_L=0.2, #Library size scale
    phi=0.1, #common dispersion
    df0=7, #BCV DoF
    x0=1, #dropout modpoint
    k_shape=0.5, #Dropout shape
    #Sample-level inputs
    m_mean=500, #Mean number of cells per subject
    #DNSMI Inputs
    alpha=0.35, 
    Rsq=0.85,
    N=50,
    rind=0.85,
    nx=15,
    nxPrime=20,
    ns=20,
    nxNoise=1500,
    ng=10,
    ngPrime=15,
    nh=15,
    ngNoise=1400,
    noiseSD=1,
    noiseMean=0
  )
```

This returns objects X (stacked data type 1), G (stacked data type 2), subjectIDs (vector specifying which rows correspond to which subjects), and y (vector of outcomes).

## STEP2: Generate similarity tensor

In order to perform the tensor regression, the similarity matrix must be estimated for each subject to be used as the predictor tensor. This will be done using the estimateSimilarity function from the estimateSimilarity.rda file.

OPTION 1: If outLoc_ was specified in STEP1 to save subject-level data in separate files:

```
  library(rTensor)
  source("estimateSimilarity.R")
  
  y <- readRDS(paste0(outputLocation, "/Ysim.rda"))
  Z_tensor <- estimateSimilarity(y=y,
                                 outLoc_=outputLocation,
                                 zeroAsMissing=FALSE)
```

OPTION 2: Otherwise, if outLoc_ was NOT specified in STEP1, then estimate the predictor tensor from the generated objects:

```
  library(rTensor)
  source("estimateSimilarity.R")
  
  y <- testNet$y
  Z_tensor <- estimateSimilarity(X=testNet$X,
                                 G=testNet$G, 
                                 y=y, 
                                 subjectIDs=testNet$subjectIDs,
                                 zeroAsMissing=FALSE)
```

Note that zeroAsMissing was set to FALSE in the above example, but if one wanted to treat zeros in the data as missing and making them be ignored in the similarity estimation, it may be updated to TRUE.

## STEP3: Tune the hyperparameters

MOSCATO applies elastic net constraints for each of the data types X and G. This requires a total of 4 hyperparameters to be tuned:

- alphaX: The weight to put on the l1-norm for X
- alphaG: The weight to put on the l1-norm for G
- maxX: The maximum number of features to be selected for X
- maxG: The maximum number of features to be selected for G

The parameters are tuned using a variation of the StARS tuning method proposed by Liu et al in 2010. numIterForStability specifies the number of subsamples to use, sizeOfSubsampleForStability specifies the number of subjects to use per subsample, and phi specifies the instability threshold. minSelect specifies the minimum number of features in X and G (vector in that order) and alphaGrid specifies the alphas to consider when tuning alphaX and alphaG.

Note that the tuning takes a long time to process, possibly days.

```
  source("tuneMOSCATO.R")
  source("supportingMOSCATOFunctions.R")  
  library(rTensor)
  library(parallel)
  
  tuned <- tuneMOSCATO(Z_tensor=Z_tensor,
                       y=y,
                       distY="gaussian",
                       numIterForStability=20,
                       sizeOfSubsampleForStability=20,
                       phi=0.01,
                       minSelect=c(5,5),
                       alphaGrid=c(0.3,0.7))
                       
  #Tuned values?
  alphaX <- tuned$tunedAlpha[1]
  alphaG <- tuned$tunedAlpha[2]
  maxX <- tuned$tunedMax[1]
  maxG <- tuned$tunedMax[2]
```

## STEP4: Perform MOSCATO with tuned values

Using the tuned alpha and max values from STEP3, the final feature selections may be made.

```
  source("supportingMOSCATOFunctions.R")    
  library(rTensor)
  
  MOSCATOFinal <- tensorElasticNet(Z_tensor=Z_tensor,
                                   y=y,
                                   glmnetAlpha=tuned$tunedAlpha,
                                   maxSelect=tuned$tunedMax
  )
  
  #What features were selected?
  selectX <- which(MOSCATOFinal$returnX==1)
  selectG <- which(MOSCATOFinal$returnG==1)
```

The tensorElasticNet function may return either "Did not converge!" or "Converged!" in the log. If "Did not converge!" was returned, this means that the Block Relaxation Algorithm did not converge on a consistent solution in estimating the coefficients from the tensor regression model. Sometimes re-running the method resolves this issue due to different random initializations.
