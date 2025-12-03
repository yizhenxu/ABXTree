#'Predict Method for BART Fits
#'
#'Predicted values based on Bayesian Additive Regression Trees model object,
#'@param obj Fitted model object from BART_mod,
#'@param newdata The data matrix to look for variables with which to predict, typically without the response; if Gcomp equals TRUE, number of rows should equal n x ndraws, where n is the number of subject. For the ith subject, (i, i+n, ..., i+(ndraws-1)n) rows are its simulated posterior outcomes from the previous simulation,
#'@param Gcomp Make predictions in the format of dynamic G-computation if TRUE; the default is FALSE,
#'@param nthin Number of posterior samples to skip between every two draws,
#'@return treefit ndraws x n posterior matrix of the sum of trees fit,
#'@return samp_y ndraws x n posterior matrix of the simulated outcome,
#'@return samp_treefit ndraws x n posterior matrix of predicted outcome based only on the sum of trees fit; this is only for type equals "multinomial",
#'@examples
#'##simulate data (example from Friedman MARS paper)
#'f = function(x){
#'  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#'}
#'sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
#'n = 100 #number of observations
#'set.seed(99)
#'x=matrix(runif(n*10),2*n,10) #10 variables, only first 5 matter
#'Ey = f(x)
#'y=Ey+sigma*rnorm(n)
#'tmp = data.frame(x,y)
#'dat = tmp[1:n,];newx = tmp[(n+1):(2*n), 1:10]
#'fml = as.formula("y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10")
#'bfit = model_bart(fml, data = dat, type = "continuous",
#'                 Prior = list(nu = 3, sigq = 0.9,
#'                              ntrees = 100,
#'                              kfac = 2,
#'                              pswap = 0.1, pbd = 0.5, pb = 0.25,
#'                              alpha = 0.95, beta = 2.0,
#'                              nc = 100, minobsnode = 10),
#'                 Mcmc = list(burn=100, ndraws = 1000))
#'
#'newy = predict_bart(obj = bfit, nthin = 0, newdata = newx)
#'@import stats
#'@export
#'@useDynLib GcompBART
predict_bart  <- function(obj, nthin = NULL, newdata = NULL, Gcomp = FALSE)
{
  if(is.data.frame(newdata)) newdata = as.matrix(newdata)
  
  if (!is.null(newdata)){
    xcolnames = obj$xcolnames
    
    if(length(xcolnames) == 1 ){
      Xtest <- data.frame(newdata[,xcolnames])
      names(Xtest) <- xcolnames[1]
    } else {
      Xtest <- newdata[,xcolnames]
    }
    
    testn = nrow(Xtest)
  } else {
    stop("newdata can not be NULL.")
  }
  
  if(is.null(nthin)){
    nthin = 0
  }
  
  ntrees = obj$ntrees
  burn = obj$burn
  npost = obj$ndraws
  
  dloc = seq(1, npost, nthin+1)
  ndraws = length(dloc)
  if(obj$type %in% c("continuous", "binary")){
    sigmasample = obj$sigmasample[dloc] #on the scale of original data
  } else {
    stop("check type")
  }
  
  if(Gcomp){
    testn = testn / ndraws
    if(testn %% 1 != 0)
      stop(paste("Error: testn does not equal to n x ndraws (after thinning)."))
  }
  
  cat("\nNumber of trees: ", paste(ntrees," "), ".\n\n", sep="")
  cat("burn-in: ", burn, ".\n\n", sep="")
  cat("Number of posteriors after burn-in: ", npost, ".\n\n", sep="")
  cat("Number of draws after thinning: ", ndraws, ".\n\n", sep="")
  cat("Number of samples: ", testn, ".\n\n", sep="")
  
  L1 = obj$TreeMod[[1]]
  L2 = obj$TreeMod[[2]]
  L3 = obj$TreeMod[[3]]
  L4 = obj$TreeMod[[4]]
  L5 = obj$TreeMod[[5]]
  L6 = obj$TreeMod[[6]]
  L7 = obj$TreeMod[[7]]
  L8 = obj$TreeMod[[8]]

  xi =  matrix(unlist(obj$xi), ncol = length(obj$xi[[1]]), byrow = TRUE) 
  
  if(obj$type == "continuous"){
    rgy = obj$rgy
    
    res =   mybartpred(as.integer(Gcomp),
                       L1,
                       L2,
                       L3,
                       L4,
                       L5,
                       L6,
                       L7,
                       L8,
                       Xtest,
                       as.integer(ncol(Xtest)),
                       as.integer(testn), as.integer(nthin),
                       as.integer(ndraws), as.integer(npost),
                       as.integer(burn),
                       as.integer(ntrees),
                       xi)
    
    #vec_test is the sum of tree fits
    
    treefit = (rgy[2]-rgy[1])*(res$vec_test+.5) + rgy[1]
    stmp = rep(sigmasample,each = testn)
    samp_y = rnorm(testn*ndraws, treefit, stmp)
    samp_y = matrix(samp_y, nrow = testn)
    treefit = matrix(treefit, nrow = testn)

  } else if(obj$type == "binary"){
    
    res =   mybartpred(as.integer(Gcomp),
                       L1,
                       L2,
                       L3,
                       L4,
                       L5,
                       L6,
                       L7,
                       L8,
                       Xtest,
                       as.integer(ncol(Xtest)),
                       as.integer(testn), as.integer(nthin),
                       as.integer(ndraws), as.integer(npost),
                       as.integer(burn),
                       as.integer(ntrees),
                       xi)
    treefit = res$vec_test
    samp_y = rnorm(testn*ndraws, treefit, 1)
    samp_y = 1*(samp_y > 0)
    samp_y = matrix(samp_y, nrow = testn)
       
    treefit = matrix(treefit, nrow = testn)
    
  } 
    
  
  ret = list(treefit = treefit,
             samp_y = samp_y);
  
  
  
  return(ret)
  
}
