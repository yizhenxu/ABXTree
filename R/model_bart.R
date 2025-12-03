#'Bayesian Additive Regression Trees
#'
#'Bayesian Additive Regression Trees Modeling for Continuous, Binary, or Multinomial Outcome,
#'@param formula response ~ covariates,
#'@param Yname Name of response,
#'@param Xname List of covariate names, with the ith element being the vector of covariate names for the ith latent variable; the order of the latent variables is in the numeric order of response categories after excluding the reference level,
#'@param data Training Data with the response and covariates in formula or Xname,
#'@param type Type of response, options are "continuous", "binary", and "multinomial",
#'@param base In the case of multinomial response, sets the reference level of the multinomial response. For example, if the response takes values 2, 3, and 4, then base = "3" sets response value 3 as the reference. Default is the highest class,
#'@param Prior List of Priors for continuous BART: e.g., Prior = list(nu=3,sigq=0.9, ntrees=200,  kfac=2.0, pswap=0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, minobsnode = 10). List of Priors for Probit BART: e.g., Prior = list(ntrees=200,  kfac=2.0, pswap=0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, minobsnode = 10).
#'The components of Prior are
#' \itemize{
#'\item nu: For continuous response: the degree of freedom in the inverse chi-square prior distribution of the error variance.
#'\item sigq: In the case of continuous response, the quantile of the error variance prior that the rough estimate (from linear regression) is placed at.
#'\item ntrees : The total number of trees in each round of BART fitting, or a vector with the ith element being the number of trees for the ith latent variable when Yname and Xname are used,
#'\item kfac : A tuning parameter that satisfies mu - kfac * sigma = ymin and mu + kfac * sigma = ymax, where mu and sigma are the mean and std of the Gaussian prior distribution of the sum of fit of all trees.
#'\item pswap : The prior probability of swap move in simulating trees; default 0, there should be pswap no larger than 1-pbd.
#'\item pbd : The prior probability of birth/death move in simulating trees; default 1.
#'\item pb : The prior probability of birth given birth/death; default 0.5.
#'\item alpha : The prior probability of a bottom node splits is alpha/(1+d)^beta, d is depth of node.
#'\item beta : see alpha.
#'\item nc : The number of equally spaced cutpoints between min and max for each covariate.
#'\item minobsnode : The minimum number of observations in bottom nodes for birth in simulating trees.
#'}
#'@param Mcmc List of MCMC starting values, burn-in ... for continuous or binary response: e.g., list(burn = 100, ndraws = 1000); for multinomial response: e.g.  list(sigma0 = diag(p - 1), burn = 100, ndraws = 1000, nSigDr = 50)
#'#'The components of Mcmc are
#' \itemize{
#'\item w0 : The n x nlatent matrix of starting values for the latent variables; this can be results from the MNP. 
#'\item sigma0 : The nlatent x nlatent matrix of starting value for the covariance matrix of latent variables.
#'\item nSigDr: User-specified upper limit to repeated draws of the covariance variance matrix of latent variables in each round of posterior draw when condition 10 in Jiao and van Dyk 2015 is not satisfied. Default 50. 
#'}
#'@param diagnostics Returns convergence diagnostics and variable inclusion proportions if True (default),
#'@param make.prediction Returns simulated outcome samp_y if TRUE (default); FALSE is only applicable to continuous and binary outcomes,
#'@return treefit ndraws x n posterior matrix of the training data sum of trees fit,
#'@return samp_y ndraws x n posterior matrix of the simulated outcome,
#'@return sigmasample posterior samples of the error standard deviation.
#'@return Percent_Acceptance Percent acceptance of Metropolis-Hastings proposals across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Nodes Average number of tree nodes across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Leaves Average number of leaves across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Depth Average tree depth across the ntrees number of trees for each posterior draw after burn-in,
#'@return Inclusion_Proportions Predictor inclusion frequencies. Smaller value of ntrees (such as 10, 20, 50, 100) is recommended for the purposes of variable selection.
#'@return TreeMod Object of the estimated sum of trees model; Treemod and sigmasample are needed for reusing the estimated model for predictions. 
#'@examples
#'##simulate data (example from Friedman MARS paper)
#'f = function(x){
#'  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#'}
#'sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
#'n = 100 #number of observations
#'set.seed(99)
#'x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
#'Ey = f(x)
#'y=Ey+sigma*rnorm(n)
#'dat = data.frame(x,y)
#'fml = y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10
#'bfit = model_bart(fml, data = dat, type = "continuous",
#'                 Prior = list(nu = 3, sigq = 0.9,
#'                              ntrees = 100,
#'                              kfac = 2,
#'                              pswap = 0.1, pbd = 0.5, pb = 0.25,
#'                              alpha = 0.95, beta = 2.0,
#'                              nc = 100, minobsnode = 10),
#'                 Mcmc = list(burn=100, ndraws = 1000))
#'
#'@import stats mlogit
#'@export
#'@useDynLib GcompBART
model_bart  <- function(formula = NULL, Yname = NULL, Xname = NULL, data, type, base = NULL,
                        Prior = NULL, Mcmc = NULL, diagnostics = TRUE, 
                        make.prediction = TRUE)
{
  
  if(!is.null(formula)){
    
    callT <- match.call(expand.dots = TRUE)
    
    formula <- mFormula(formula) # mFormula is from package mlogit
    
    response.name <- paste(deparse(attr(formula, "lhs")[[1L]]))
    
    # mf = list(function name, formula, data) 
    m <- match(c("formula", "data"), names(callT), 0L)
    mf <- callT
    mf <- mf[c(1L, m)]
    
    mf$formula <- formula
    
    # mf = model.frame(formula, data), 
    # returns data frame of attr(terms(mf),'variables'), which is (y, covariates)
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    # extract response from the mf data frame based on attr(terms(mf),'response'), col no. of y (1)
    Y <- model.response(mf)
    
    Terms <- attr(mf, "terms") # <-> terms(mf)
    
    # extract design matrix from mf based on Terms, which has attr(terms(mf),'intercept') = 1
    X <- model.matrix.default(Terms, mf) 
    xcolnames <- colnames(X)[-1]
    # this is equiv to doing:
    # attr(Terms,'intercept') = 0; X = ...; xcolnames = colnames(X)
    # X <- X[,xcolnames] would be unnecessary
    
    if(length(xcolnames) == 1 ){
      X <- data.frame(X[,xcolnames])
      names(X) <- xcolnames[1]
    } else {
      X <- X[,xcolnames]
    }
  } else {
    
    ncov = c()
    data = as.data.frame(data)
    
    lXname = length(Xname)
    fmllist = vector("list",lXname)
    xcolnamesL = vector("list",lXname)
    for(j in 1:lXname){
      
      # for each latent, generate the corresponding formula
      Xj = Xname[[j]]
      lX = length(Xj)
      if(lX == 0) stop("Error: Each element of X should include at least one covariate.")
      if(lX == 1){
        fml = paste0(Yname,"~", Xj[lX],collapse = "")
      } else {
        fml = paste0(Yname,"~",paste0(Xj[-lX],"+",collapse = ""), Xj[lX],collapse = "")
      }
      fmllist[[j]] = as.formula(fml)
      
      # extract design matrix
      mf = model.frame(fml, data)
      Terms = attr(mf, "terms")
      if(j==1){
        X = model.matrix.default(Terms, mf)
        xcolnames = xcolnamesL[[j]] = colnames(X)[-1]  
        ncov = c(ncov, length(xcolnames))
      } else {
        Xnew = model.matrix.default(Terms, mf)
        X = cbind(X, Xnew)
        xcolnamesL[[j]] = colnames(Xnew)[-1]
        xcolnames = c(xcolnames, xcolnamesL[[j]])
        ncov = c(ncov, length(colnames(Xnew)[-1]))
      }
      
    }
    
    X = X[, xcolnames]
    
    Y <- model.response(mf)
    
  }
  
  n = nrow(X)
  
  binaryX = rep(NA,ncol(X))
  for(i in 1:ncol(X)){
    binaryX[i] = 1*(length(unique(X[,i]))==2)
  }
  
  if(missing(Prior))
  {
    ntrees=200; kfac=2.0;pswap=0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; minobsnode = 10;
  }
  else
  {
    if(is.null(Prior$ntrees)) {ntrees=200} else {ntrees=Prior$ntrees}
    if(is.null(Prior$kfac)) {kfac=2.0} else {kfac=Prior$kfac}
    if(is.null(Prior$pswap)) {pswap=0.0} else {pswap=Prior$pswap}
    if(is.null(Prior$pbd)) {pbd=1.0} else {pbd=Prior$pbd}
    if(is.null(Prior$pb)) {pb=0.5} else {pb=Prior$pb}
    if(is.null(Prior$beta)) {beta = 2.0} else {beta=Prior$beta}
    if(is.null(Prior$alpha)) {alpha = 0.95} else {alpha=Prior$alpha}
    if(is.null(Prior$nc)) {nc=100} else {nc=Prior$nc}
    if(is.null(Prior$minobsnode)) {minobsnode= 10} else {minobsnode=Prior$minobsnode}
  }
  
  if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
  if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}
  
  if(type == "continuous"){
    rgy = range(Y)
    Y = -.5 + (Y-rgy[1])/(rgy[2]-rgy[1])
    Data = list(y = Y,X = X)
    
    ###
    if(missing(Prior))
    {
      nu=3;sigq=0.90;
    }
    else
    { 
      if(is.null(Prior$nu)) {nu = 3} else {nu=Prior$nu}
      if(is.null(Prior$sigq)) {sigq=0.90} else {sigq=Prior$sigq}
      
    }
    
    ###
    templm = lm(Y~X)
    sigest = summary(templm)$sigma# square root of the estimatedvariance of the random error
    sigmasample = as.double(rep(sigest, ndraws));
    
    # Pr(sigma < k) = sigq <=>
    # lambda = qchisq(1-sigq,nu)*k^2/nu
    qchi = qchisq(1.0-sigq, nu);
    lambda = (sigest*sigest*qchi)/nu;
    ###
    
  } else if(type == "binary"){
    Data = list(y = Y,X = X)
    
  } 
  
  cat("Number of trees: ", paste(ntrees," "), ".\n\n", sep="")
  cat("Number of draws: ", ndraws, ".\n\n", sep="")
  cat("burn-in: ", burn, ".\n\n", sep="") 
  
  if(type == "continuous"){
    
    res =   mybartmod(Data$X, 
                      as.double(sigest),
                      as.integer(n),
                      Data$y,
                      as.integer(ncol(Data$X)),
                      as.integer(nu),
                      as.double(lambda),
                      as.integer(ndraws),
                      as.integer(burn),
                      as.integer(ntrees),
                      as.double(kfac),
                      as.double(pswap),
                      as.double(pbd),
                      as.double(pb),
                      as.double(alpha),
                      as.double(beta),
                      as.integer(nc),
                      as.integer(minobsnode),
                      binaryX,
                      as.integer(diagnostics))
    
    names(res$Inclusion_Proportions) = xcolnames
    
    res$sigmasample = res$sigmasample*(rgy[2]-rgy[1])
    res$treefit = (rgy[2]-rgy[1])*(res$treefit+.5) + rgy[1] #nsamp X npost
    
    if(make.prediction){
      stmp = rep(res$sigmasample,each = n)
      res$samp_y = rnorm(n*ndraws, res$treefit, stmp)
      res$samp_y = matrix(res$samp_y, nrow = n)
    }
    
    res$treefit = matrix(res$treefit, nrow = n)
    
    res = append(res, list("xcolnames"=xcolnames,
                           "rgy" = rgy,
                           "ntrees" = ntrees,
                           "burn" = burn,
                           "ndraws" = ndraws,
                           "type" = "continuous"), length(res))
    
  } else if(type == "binary"){
    res =   mypbartmod(Data$X, 
                       as.integer(n),
                       Data$y,
                       as.integer(ncol(Data$X)),
                       as.integer(ndraws),
                       as.integer(burn),
                       as.integer(ntrees),
                       as.double(kfac),
                       as.double(pswap),
                       as.double(pbd),
                       as.double(pb),
                       as.double(alpha),
                       as.double(beta),
                       as.integer(nc),
                       as.integer(minobsnode),
                       binaryX,
                       as.integer(diagnostics))
    
    names(res$Inclusion_Proportions) = xcolnames
    
    if(make.prediction){
      res$samp_y = rnorm(n*ndraws, res$treefit, 1)
      res$samp_y = 1*(res$samp_y > 0)
      res$samp_y = matrix(res$samp_y, nrow = n)
    }
    
    # sum of tree fits
    res$treefit = matrix(res$treefit, nrow = n)
    
    res = append(res, list("xcolnames"=xcolnames,
                           "ntrees" = ntrees,
                           "burn" = burn,
                           "ndraws" = ndraws,
                           "type" = "binary"), length(res))
    
  }    
  
  
  
  return(res)
  
  
}
