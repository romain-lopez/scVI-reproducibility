optimleft_fun <- function(i, gamma_mu, gamma_pi, W, Y, V_mu, alpha_mu,
                          X_mu, beta_mu, O_mu, V_pi, alpha_pi, X_pi,
                          beta_pi, O_pi, zeta, epsilonleft) {
  optim( fn=zinb.loglik.regression,
         gr=zinb.loglik.regression.gradient,
         par=c(gamma_mu[,i], gamma_pi[,i],
               t(W[i,])),
         Y=t(Y[i,]),
         A.mu=V_mu,
         B.mu=t(alpha_mu),
         C.mu=t(X_mu[i,]%*%beta_mu + O_mu[i,]),
         A.pi=V_pi,
         B.pi=t(alpha_pi),
         C.pi=t(X_pi[i,]%*%beta_pi + O_pi[i,]),
         C.theta=zeta,
         epsilon=epsilonleft,
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}

optimright_fun <- function(j, beta_mu, alpha_mu, beta_pi, alpha_pi,
                           Y, X_mu, W, V_mu, gamma_mu, O_mu, X_pi,
                           V_pi, gamma_pi, O_pi, zeta, n, epsilonright) {
  optim( fn=zinb.loglik.regression,
         gr=zinb.loglik.regression.gradient,
         par=c(beta_mu[,j], alpha_mu[,j],
               beta_pi[,j], alpha_pi[,j]),
         Y=Y[,j], A.mu=cbind(X_mu, W),
         C.mu=t(V_mu[j,] %*% gamma_mu) + O_mu[,j],
         A.pi=cbind(X_pi, W),
         C.pi=t(V_pi[j,] %*% gamma_pi) + O_pi[,j],
         C.theta=matrix(zeta[j], nrow = n, ncol = 1),
         epsilon=epsilonright,
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}



#' TRICKILY Optimize the cell parameters of a ZINB regression model
#'
#' The parameters of the model W and gamma given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument. It is recommended
#' to call zinb_initialize and zinb_optimize before calling this function as this will only
#' only optimize for a cross-validation set. 
#' @param m_train The fitted model of class ZinbModel on train set
#' @param m_val The init model of class ZinbModel on val set
#' @param Y_val The validation matrix of counts.
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters gamma_mu, gamma_pi and W fitted on a validation set.
#' @examples
#' Y_train = matrix(10, 3, 5)
#' Y_val = matrix(10, 3, 5)
#' m_train = zinbModel(n=NROW(Y_train), J=NCOL(Y_train))
#' m_train = zinbInitialize(m_train, Y_train)
#' m_val = zinbModel(n=NROW(Y_val), J=NCOL(Y_val))
#' m_val = zinbInitialize(m_val, Y_val)
#
#' m_train = zinbOptimize(m, Y_train)
#' m_val = zinbOptimizeVal(m_train, m_val, Y_val)
#' 
#' @export
zinbOptimizeValTrick <- function(m_train, m_val, Y_val, maxiter=25,
                            stop.epsilon=.0001, verbose=FALSE,
                            BPPARAM=BiocParallel::bpparam()) {
  
  #' FIRST STEP: copy shared parameters from m_train and static param from 
  #' m_val to create m to be optimized
  
  m <- zinbModel(
    # Fixed values shared by models
    X = m_val@X, V = m_val@V, O_mu = m_val@O_mu,
    O_pi = m_val@O_pi, which_X_mu = m_val@which_X_mu,
    which_X_pi = m_val@which_X_pi, which_V_mu = m_val@which_V_mu,
    which_V_pi = m_val@which_V_pi, 
    
    # Values initialized with a different shape
    W = m_val@W, 
    gamma_mu = m_val@gamma_mu, gamma_pi = m_val@gamma_pi,
    
    # Values from training set
    beta_mu = m_train@beta_mu, beta_pi = m_train@beta_pi,
    alpha_mu = m_train@alpha_mu, alpha_pi = m_train@alpha_pi, 
    zeta = m_train@zeta, 
    
    # Fixed values also shared
    epsilon_beta_mu = m_train@epsilon_beta_mu,
    epsilon_gamma_mu = m_train@epsilon_gamma_mu,
    epsilon_beta_pi = m_train@epsilon_beta_pi,
    epsilon_gamma_pi = m_train@epsilon_gamma_pi,
    epsilon_W = m_train@epsilon_W, 
    epsilon_alpha = m_train@epsilon_alpha,
    epsilon_zeta = m_train@epsilon_zeta,
    epsilon_min_logit = m_train@epsilon_min_logit)
  
  #' now optimize the RIDGE from zinbOptimize
  total.lik=rep(NA,maxiter)
  n <- nSamples(m)
  J <- nFeatures(m)
  
  
  epsilonleft <- c(getEpsilon_gamma_mu(m),
                   getEpsilon_gamma_pi(m), getEpsilon_W(m))
  nleft <- c(length(getEpsilon_gamma_mu(m)),
             length(getEpsilon_gamma_pi(m)), length(getEpsilon_W(m)))
  optimleft = (sum(nleft)>0)
  
  # extract fixed quantities from m
  X_mu <- getX_mu(m)
  V_mu <- getV_mu(m)
  X_pi <- getX_pi(m)
  V_pi <- getV_pi(m)
  O_mu <- m@O_mu
  O_pi <- m@O_pi
  
  # exctract paramters from m (remember to update!)
  beta_mu <- getBeta_mu(m)
  alpha_mu <- getAlpha_mu(m)
  gamma_mu <- 0 * getGamma_mu(m)
  gamma_pi <- 0 * getGamma_pi(m)
  beta_pi <- getBeta_pi(m)
  alpha_pi <- getAlpha_pi(m)
  W <- getW(m)
  zeta <- getZeta(m)
  
  for (iter in seq_len(maxiter)){
    if (verbose) {message("Iteration ",iter)}
    
    # Evaluate total penalized likelihood
    mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                W %*% alpha_mu + O_mu)
    
    logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
      W %*% alpha_pi + O_pi
    
    theta <- exp(zeta)
    
    loglik <- zinb.loglik(Y_val, mu, rep(theta, rep(n, J)), logitPi)
    
    penalty <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
      sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
      sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
      sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
      sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
      sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
      sum(getEpsilon_W(m)*t(W)^2)/2 +
      getEpsilon_zeta(m)*var(zeta)/2
    
    total.lik[iter] <- loglik - penalty
    
    if (verbose) {message("penalized log-likelihood = ",
                          total.lik[iter])}
    
    # If the increase in likelihood is smaller than 0.5%, stop maximization
    if(iter > 1){
      if(abs((total.lik[iter]-total.lik[iter-1]) /
             total.lik[iter-1])<stop.epsilon)
        break
    }
    
    # Optimize left factors
    if (optimleft) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(n), optimleft_fun,
                 gamma_mu, gamma_pi, W, Y_val, V_mu, alpha_mu,
                 X_mu, beta_mu, O_mu, V_pi, alpha_pi, X_pi,
                 beta_pi, O_pi, zeta, epsilonleft,
                 BPPARAM=BPPARAM)), nrow=sum(nleft))
      
      if (verbose) {print(proc.time()-ptm)}
      ind <- 1
      if (nleft[1]>0) {
        #gamma_mu <- estimate[ind:(ind+nleft[1]-1),,drop=FALSE]
        ind <- ind+nleft[1]
      }
      if (nleft[2]>0) {
        #gamma_pi <- estimate[ind:(ind+nleft[2]-1),,drop=FALSE]
        ind <- ind+nleft[2]
      }
      if (nleft[3]>0) {
        W <- t(estimate[ind:(ind+nleft[3]-1),,drop=FALSE])
        ind <- ind+nleft[3]
      }
    }
  }
  
  out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu, O_pi = m@O_pi,
                   which_X_mu = m@which_X_mu, which_X_pi = m@which_X_pi,
                   which_V_mu = m@which_V_mu, which_V_pi = m@which_V_pi,
                   W = W, beta_mu = beta_mu, beta_pi = beta_pi,
                   gamma_mu = gamma_mu, gamma_pi = gamma_pi,
                   alpha_mu = alpha_mu, alpha_pi = alpha_pi, zeta = zeta,
                   epsilon_beta_mu = m@epsilon_beta_mu,
                   epsilon_gamma_mu = m@epsilon_gamma_mu,
                   epsilon_beta_pi = m@epsilon_beta_pi,
                   epsilon_gamma_pi = m@epsilon_gamma_pi,
                   epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                   epsilon_zeta = m@epsilon_zeta,
                   epsilon_min_logit = m@epsilon_min_logit)
  
  return(out)
}


#' Optimize the cell parameters of a ZINB regression model
#'
#' The parameters of the model W and gamma given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument. It is recommended
#' to call zinb_initialize and zinb_optimize before calling this function as this will only
#' only optimize for a cross-validation set. 
#' @param m_train The fitted model of class ZinbModel on train set
#' @param m_val The init model of class ZinbModel on val set
#' @param Y_val The validation matrix of counts.
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters gamma_mu, gamma_pi and W fitted on a validation set.
#' @examples
#' Y_train = matrix(10, 3, 5)
#' Y_val = matrix(10, 3, 5)
#' m_train = zinbModel(n=NROW(Y_train), J=NCOL(Y_train))
#' m_train = zinbInitialize(m_train, Y_train)
#' m_val = zinbModel(n=NROW(Y_val), J=NCOL(Y_val))
#' m_val = zinbInitialize(m_val, Y_val)
#
#' m_train = zinbOptimize(m, Y_train)
#' m_val = zinbOptimizeVal(m_train, m_val, Y_val)
#' 
#' @export
zinbOptimizeVal <- function(m_train, m_val, Y_val, maxiter=25,
                         stop.epsilon=.0001, verbose=FALSE,
                         BPPARAM=BiocParallel::bpparam()) {
  
  #' FIRST STEP: copy shared parameters from m_train and static param from 
  #' m_val to create m to be optimized
  
  m <- zinbModel(
                   # Fixed values shared by models
                   X = m_val@X, V = m_val@V, O_mu = m_val@O_mu,
                   O_pi = m_val@O_pi, which_X_mu = m_val@which_X_mu,
                   which_X_pi = m_val@which_X_pi, which_V_mu = m_val@which_V_mu,
                   which_V_pi = m_val@which_V_pi, 
                   
                   # Values initialized with a different shape
                   W = m_val@W, 
                   gamma_mu = m_val@gamma_mu, gamma_pi = m_val@gamma_pi,
                   
                   # Values from training set
                   beta_mu = m_train@beta_mu, beta_pi = m_train@beta_pi,
                   alpha_mu = m_train@alpha_mu, alpha_pi = m_train@alpha_pi, 
                   zeta = m_train@zeta, 
                   
                   # Fixed values also shared
                   epsilon_beta_mu = m_train@epsilon_beta_mu,
                   epsilon_gamma_mu = m_train@epsilon_gamma_mu,
                   epsilon_beta_pi = m_train@epsilon_beta_pi,
                   epsilon_gamma_pi = m_train@epsilon_gamma_pi,
                   epsilon_W = m_train@epsilon_W, 
                   epsilon_alpha = m_train@epsilon_alpha,
                   epsilon_zeta = m_train@epsilon_zeta,
                   epsilon_min_logit = m_train@epsilon_min_logit)
  
  #' now optimize the RIDGE from zinbOptimize
  total.lik=rep(NA,maxiter)
  n <- nSamples(m)
  J <- nFeatures(m)

  
  epsilonleft <- c(getEpsilon_gamma_mu(m),
                   getEpsilon_gamma_pi(m), getEpsilon_W(m))
  nleft <- c(length(getEpsilon_gamma_mu(m)),
             length(getEpsilon_gamma_pi(m)), length(getEpsilon_W(m)))
  optimleft = (sum(nleft)>0)
  
  # extract fixed quantities from m
  X_mu <- getX_mu(m)
  V_mu <- getV_mu(m)
  X_pi <- getX_pi(m)
  V_pi <- getV_pi(m)
  O_mu <- m@O_mu
  O_pi <- m@O_pi
  
  # exctract paramters from m (remember to update!)
  beta_mu <- getBeta_mu(m)
  alpha_mu <- getAlpha_mu(m)
  gamma_mu <- getGamma_mu(m)
  gamma_pi <- getGamma_pi(m)
  beta_pi <- getBeta_pi(m)
  alpha_pi <- getAlpha_pi(m)
  W <- getW(m)
  zeta <- getZeta(m)
  
  for (iter in seq_len(maxiter)){
    if (verbose) {message("Iteration ",iter)}
    
    # Evaluate total penalized likelihood
    mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                W %*% alpha_mu + O_mu)
    
    logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
      W %*% alpha_pi + O_pi
    
    theta <- exp(zeta)
    
    loglik <- zinb.loglik(Y_val, mu, rep(theta, rep(n, J)), logitPi)
    
    penalty <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
      sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
      sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
      sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
      sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
      sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
      sum(getEpsilon_W(m)*t(W)^2)/2 +
      getEpsilon_zeta(m)*var(zeta)/2
    
    total.lik[iter] <- loglik - penalty
    
    if (verbose) {message("penalized log-likelihood = ",
                          total.lik[iter])}
    
    # If the increase in likelihood is smaller than 0.5%, stop maximization
    if(iter > 1){
      if(abs((total.lik[iter]-total.lik[iter-1]) /
             total.lik[iter-1])<stop.epsilon)
        break
    }
  
    # Optimize left factors
    if (optimleft) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(n), optimleft_fun,
                 gamma_mu, gamma_pi, W, Y_val, V_mu, alpha_mu,
                 X_mu, beta_mu, O_mu, V_pi, alpha_pi, X_pi,
                 beta_pi, O_pi, zeta, epsilonleft,
                 BPPARAM=BPPARAM)), nrow=sum(nleft))
      
      if (verbose) {print(proc.time()-ptm)}
      ind <- 1
      if (nleft[1]>0) {
        gamma_mu <- estimate[ind:(ind+nleft[1]-1),,drop=FALSE]
        ind <- ind+nleft[1]
      }
      if (nleft[2]>0) {
        gamma_pi <- estimate[ind:(ind+nleft[2]-1),,drop=FALSE]
        ind <- ind+nleft[2]
      }
      if (nleft[3]>0) {
        W <- t(estimate[ind:(ind+nleft[3]-1),,drop=FALSE])
        ind <- ind+nleft[3]
      }
    }
  }
  
  out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu, O_pi = m@O_pi,
                   which_X_mu = m@which_X_mu, which_X_pi = m@which_X_pi,
                   which_V_mu = m@which_V_mu, which_V_pi = m@which_V_pi,
                   W = W, beta_mu = beta_mu, beta_pi = beta_pi,
                   gamma_mu = gamma_mu, gamma_pi = gamma_pi,
                   alpha_mu = alpha_mu, alpha_pi = alpha_pi, zeta = zeta,
                   epsilon_beta_mu = m@epsilon_beta_mu,
                   epsilon_gamma_mu = m@epsilon_gamma_mu,
                   epsilon_beta_pi = m@epsilon_beta_pi,
                   epsilon_gamma_pi = m@epsilon_gamma_pi,
                   epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                   epsilon_zeta = m@epsilon_zeta,
                   epsilon_min_logit = m@epsilon_min_logit)

  return(out)
}

#' Computes likelihood as a score
#'
#' The parameters of the model W and gamma given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument. It is recommended
#' to call zinb_initialize and zinb_optimize before calling this function as this will only
#' only optimize for a cross-validation set. 
#' @param m The fitted model of class ZinbModel
#' @param Y The  matrix of counts.
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters gamma_mu, gamma_pi and W fitted on a validation set.
#' @examples
#' Y_train = matrix(10, 3, 5)
#' Y_val = matrix(10, 3, 5)
#' m_train = zinbModel(n=NROW(Y_train), J=NCOL(Y_train))
#' m_train = zinbInitialize(m_train, Y_train)
#' m_val = zinbModel(n=NROW(Y_val), J=NCOL(Y_val))
#' m_val = zinbInitialize(m_val, Y_val)
#
#' m_train = zinbOptimize(m, Y_train)
#' m_val = zinbOptimizeVal(m_train, m_val, Y_val)
#' score = zinbLikelihood(m_val, Y_val)
#' 
#' @export
zinbLikelihood <- function(m, Y, maxiter=25,
                            stop.epsilon=.0001, verbose=FALSE,
                            BPPARAM=BiocParallel::bpparam()) {
  n <- nSamples(m)
  J <- nFeatures(m)
# extract fixed quantities from m
X_mu <- getX_mu(m)
V_mu <- getV_mu(m)
X_pi <- getX_pi(m)
V_pi <- getV_pi(m)
O_mu <- m@O_mu
O_pi <- m@O_pi

# exctract paramters from m (remember to update!)
beta_mu <- getBeta_mu(m)
alpha_mu <- getAlpha_mu(m)
gamma_mu <- getGamma_mu(m)
beta_pi <- getBeta_pi(m)
alpha_pi <- getAlpha_pi(m)
gamma_pi <- getGamma_pi(m)
W <- getW(m)
zeta <- getZeta(m)

mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
            W %*% alpha_mu + O_mu)

logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
  W %*% alpha_pi + O_pi

theta <- exp(zeta)

loglik <- zinb.loglik(Y, mu, rep(theta, rep(n, J)), logitPi)
}

#' Get infered ZINB Fit once model has been fitted on the data using ZinbOptimize or ZinbOptimizeVal
#' @param m The fitted model of class ZinbModel
#' @return estimated mu on the fitted model
#' @examples
#' Y = matrix(10, 3, 5)
#' m = zinbModel(n=NROW(Y), J=NCOL(Y))
#' m = zinbInitialize(m, Y)
#' m = zinbOptimize(m, Y)
#' fit = zinbFit(m)
#' 
#' @export
zinbFit <- function(m) {
  n <- nSamples(m)
  J <- nFeatures(m)
  # extract fixed quantities from m
  X_mu <- getX_mu(m)
  V_mu <- getV_mu(m)
  X_pi <- getX_pi(m)
  V_pi <- getV_pi(m)
  O_mu <- m@O_mu
  O_pi <- m@O_pi
  
  # exctract paramters from m (remember to update!)
  beta_mu <- getBeta_mu(m)
  alpha_mu <- getAlpha_mu(m)
  gamma_mu <- getGamma_mu(m)
  beta_pi <- getBeta_pi(m)
  alpha_pi <- getAlpha_pi(m)
  gamma_pi <- getGamma_pi(m)
  W <- getW(m)
  zeta <- getZeta(m)
  
  mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
              W %*% alpha_mu + O_mu)
  
  logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
    W %*% alpha_pi + O_pi
  
  theta <- exp(zeta)
  
  res <- list(mu, logitPi, theta)
}


#' TRICKILY Optimize the parameters of a ZINB regression model
#'
#' The parameters of the model given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument. It is recommended
#' to call zinb_initialize before this function to have good starting point for
#' optimization, since the optimization problem is not convex and can only
#' converge to a local minimum.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @param commondispersion Whether the dispersion is the same for all features
#'   (default=TRUE)
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu,
#'   gamma_pi, W.
#' @examples
#' Y = matrix(10, 3, 5)
#' m = zinbModel(n=NROW(Y), J=NCOL(Y))
#' m = zinbInitialize(m, Y)
#' m = zinbOptimize(m, Y)
#' @export
zinbOptimizeTrick <- function(m, Y, commondispersion=TRUE, maxiter=25,
                         stop.epsilon=.0001, verbose=FALSE,
                         BPPARAM=BiocParallel::bpparam()) {
  
  total.lik=rep(NA,maxiter)
  n <- nSamples(m)
  J <- nFeatures(m)
  
  epsilonright <- c(getEpsilon_beta_mu(m), getEpsilon_alpha(m),
                    getEpsilon_beta_pi(m), getEpsilon_alpha(m))
  nright <- c(length(getEpsilon_beta_mu(m)), length(getEpsilon_alpha(m)),
              length(getEpsilon_beta_pi(m)), length(getEpsilon_alpha(m)))
  optimright = (sum(nright)>0)
  
  epsilonleft <- c(getEpsilon_gamma_mu(m),
                   getEpsilon_gamma_pi(m), getEpsilon_W(m))
  nleft <- c(length(getEpsilon_gamma_mu(m)),
             length(getEpsilon_gamma_pi(m)), length(getEpsilon_W(m)))
  optimleft = (sum(nleft)>0)
  
  orthog <- (nFactors(m)>0)
  
  # extract fixed quantities from m
  X_mu <- getX_mu(m)
  V_mu <- getV_mu(m)
  X_pi <- getX_pi(m)
  V_pi <- getV_pi(m)
  O_mu <- m@O_mu
  O_pi <- m@O_pi
  
  # exctract paramters from m (remember to update!)
  beta_mu <- getBeta_mu(m)
  alpha_mu <- getAlpha_mu(m)
  
  # TRICK !
  gamma_mu <- 0 * getGamma_mu(m)
  gamma_pi <- 0 * getGamma_pi(m)
  
  
  beta_pi <- getBeta_pi(m)
  alpha_pi <- getAlpha_pi(m)
  W <- getW(m)
  zeta <- getZeta(m)
  
  for (iter in seq_len(maxiter)){
    if (verbose) {message("Iteration ",iter)}
    
    # Evaluate total penalized likelihood
    mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                W %*% alpha_mu + O_mu)
    
    logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
      W %*% alpha_pi + O_pi
    
    theta <- exp(zeta)
    
    loglik <- zinb.loglik(Y, mu, rep(theta, rep(n, J)), logitPi)
    
    penalty <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
      sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
      sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
      sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
      sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
      sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
      sum(getEpsilon_W(m)*t(W)^2)/2 +
      getEpsilon_zeta(m)*var(zeta)/2
    
    total.lik[iter] <- loglik - penalty
    
    if (verbose) {message("penalized log-likelihood = ",
                          total.lik[iter])}
    
    # If the increase in likelihood is smaller than 0.5%, stop maximization
    if(iter > 1){
      if(abs((total.lik[iter]-total.lik[iter-1]) /
             total.lik[iter-1])<stop.epsilon)
        break
    }
    
    # 1. Optimize dispersion
    zeta <- zinbOptimizeDispersion(J, mu, logitPi, getEpsilon_zeta(m), Y,
                                   commondispersion=commondispersion,
                                   BPPARAM=BPPARAM)
    
    # Evaluate total penalized likelihood
    if (verbose) {
      message("After dispersion optimization = ",
              loglik(m, Y) - penalty(m))
    }
    
    # 2. Optimize right factors
    
    if (optimright) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(J), optimright_fun,
                 beta_mu, alpha_mu, beta_pi, alpha_pi,
                 Y, X_mu, W, V_mu, gamma_mu, O_mu, X_pi,
                 V_pi, gamma_pi, O_pi, zeta, n, epsilonright,
                 BPPARAM=BPPARAM)), nrow=sum(nright))
      
      if (verbose) {print(proc.time()-ptm)}
      ind <- 1
      if (nright[1]>0) {
        beta_mu <- estimate[ind:(ind+nright[1]-1),,drop=FALSE]
        ind <- ind+nright[1]
      }
      if (nright[2]>0) {
        alpha_mu <- estimate[ind:(ind+nright[2]-1),,drop=FALSE]
        ind <- ind+nright[2]
      }
      if (nright[3]>0) {
        beta_pi <- estimate[ind:(ind+nright[3]-1),,drop=FALSE]
        ind <- ind+nright[3]
      }
      if (nright[4]>0) {
        alpha_pi <- estimate[ind:(ind+nright[4]-1),,drop=FALSE]
      }
    }
    # Evaluate total penalized likelihood
    if (verbose) {
      message("After right optimization = ",
              loglik(m, Y)-penalty(m))
    }
    
    # 3. Orthogonalize
    if (orthog) {
      o <- orthogonalizeTraceNorm(W, cbind(alpha_mu,
                                           alpha_pi),
                                  m@epsilon_W, m@epsilon_alpha)
      W <- o$U
      alpha_mu <- o$V[,1:J,drop=FALSE]
      alpha_pi <- o$V[,(J+1):(2*J),drop=FALSE]
    }
    
    # Evaluate total penalized likelihood
    if (verbose) {message("After orthogonalization = ",
                          loglik(m, Y) - penalty(m))}
    
    # 4. Optimize left factors
    if (optimleft) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(n), optimleft_fun,
                 gamma_mu, gamma_pi, W, Y, V_mu, alpha_mu,
                 X_mu, beta_mu, O_mu, V_pi, alpha_pi, X_pi,
                 beta_pi, O_pi, zeta, epsilonleft,
                 BPPARAM=BPPARAM)), nrow=sum(nleft))
      
      if (verbose) {print(proc.time()-ptm)}
      ind <- 1
      if (nleft[1]>0) {
        #TRICK
        #gamma_mu <- estimate[ind:(ind+nleft[1]-1),,drop=FALSE]
        ind <- ind+nleft[1]
      }
      if (nleft[2]>0) {
        # TRICK
        #gamma_pi <- estimate[ind:(ind+nleft[2]-1),,drop=FALSE]
        ind <- ind+nleft[2]
      }
      if (nleft[3]>0) {
        W <- t(estimate[ind:(ind+nleft[3]-1),,drop=FALSE])
        ind <- ind+nleft[3]
      }
    }
    
    # Evaluate total penalized likelihood
    if (verbose) {message("After left optimization = ",
                          loglik(m, Y) - penalty(m))}
    
    # 5. Orthogonalize
    if (orthog) {
      o <- orthogonalizeTraceNorm(W, cbind(alpha_mu,
                                           alpha_pi),
                                  m@epsilon_W, m@epsilon_alpha)
      W <- o$U
      alpha_mu <- o$V[,1:J,drop=FALSE]
      alpha_pi <- o$V[,(J+1):(2*J),drop=FALSE]
    }
    # Evaluate total penalized likelihood
    if (verbose) {message("After orthogonalization = ",
                          loglik(m, Y) - penalty(m))}
    
  }
  
  out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu, O_pi = m@O_pi,
                   which_X_mu = m@which_X_mu, which_X_pi = m@which_X_pi,
                   which_V_mu = m@which_V_mu, which_V_pi = m@which_V_pi,
                   W = W, beta_mu = beta_mu, beta_pi = beta_pi,
                   gamma_mu = gamma_mu, gamma_pi = gamma_pi,
                   alpha_mu = alpha_mu, alpha_pi = alpha_pi, zeta = zeta,
                   epsilon_beta_mu = m@epsilon_beta_mu,
                   epsilon_gamma_mu = m@epsilon_gamma_mu,
                   epsilon_beta_pi = m@epsilon_beta_pi,
                   epsilon_gamma_pi = m@epsilon_gamma_pi,
                   epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                   epsilon_zeta = m@epsilon_zeta,
                   epsilon_min_logit = m@epsilon_min_logit)
  
  return(out)
}


