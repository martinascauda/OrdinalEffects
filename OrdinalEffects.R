library(mvtnorm)
library(MASS)
library(truncnorm)
library(OwenQ)
library(BiDAG)

# Major file containing the OSEM algorithm
source("./OSEMSource/R/ordinalScore.R")




##' rmvDAG2var(N, DAG, V):
##' a function that does the same thing as the pcalg::rmvDAG function (adapted from rmvDAG2 Luo et. all)
##' but the input DAG is not necessarily topologically ordered and the error variance is not necessarily the identity
##' @param N: number of samples to be drawn
##' @param V: error variance matrix
##' @param DAG: a graphNEL object 
##' @return a Gaussian dataset
rmvDAG2var <- function(N, DAG,V) {
  AM <- as(DAG, "matrix")
  sorted_ind <- ggm::topOrder((AM != 0))
  n <- nrow(AM)
  data <- matrix(nrow = N,ncol = n)
  for (j in sorted_ind) {
    parentnodes <- which(AM[,j] != 0)
    lp <- length(parentnodes)
    switch (as.character(lp),
            "0" = {data[,j] <- rnorm(N, mean = 0, sd=sqrt(V[j,j]))},
            "1" = {data[,j] <- rnorm(N, mean = data[,parentnodes] * AM[parentnodes,j], sd = sqrt(V[j,j]))},
            {data[,j] <- rnorm(N, mean = data[,parentnodes] %*% AM[parentnodes,j], sd = sqrt(V[j,j]))}
    )
  }
  return(data)
}

##' convertcutsToOrdinal(scaled_data, cuts):
##' a function that converts standardized Gaussian data into ordinal data according to cut points
##' @param scaled_data: Gaussian dataset with each dimension standardized
##' @param cuts: thresholds of ordinal levels
##' @return an ordinal dataset
convertcutsToOrdinal <- function(scaled_data, cuts) {
  n <- ncol(scaled_data)
  ordinal_data<-scaled_data
  for (i in c(1:n)){
    cuts_i <- cuts[i]
    temp <- cut(scaled_data[,i], simplify2array(cuts_i), labels = FALSE) - 1
    ordinal_data[,i] <- temp
  }
  colnames(ordinal_data) <- c(1:n)
  return(ordinal_data)
}


##' generateTRUEOrdinal(N, trueDAG, cuts, V):
##' a function that generate ordinal from a latent Gaussian DAG with mean 0 and error variance matrix V given vector of thresholds
##' @param N: number of samples to be drawn
##' @param V: error variance matrix
##' @param trueDAG: a graphNEL object 
##' @param cuts: thresholds of ordinal levels
##' @return an ordinal dataset
generateTRUEOrdinal <- function(N, trueDAG, cuts,V) {
  hidden_data <- rmvDAG2var(N, trueDAG,V)
  scaled_data <- scale(hidden_data, center =TRUE, scale = TRUE)
  
  # Convert the Gaussian dataset into an ordinal dataset
  ordinal_data <- convertcutsToOrdinal(scaled_data, cuts = cuts)
  
  return(ordinal_data)
  
}

##' extractCuts(scaled_data, exp_levels,concent_param): 
##' a function that that converts standardized Gaussian data into ordinal data and saves the cuts adopted ((adapted from convertToOrdinal Luo et. all)
##' @param scaled_data: Gaussian dataset with each dimension standardized
##' @param exp_levels: expected number of ordinal levels
##' @param concent_param: Dirichlet concentration parameter
##' @return an ordinal dataset and the cuts adopted to generate it 
extractCuts <- function(scaled_data, exp_levels = 4,concent_param = 2) {
  n <- ncol(scaled_data)
  if (exp_levels == 2) {
    ordinal_levels <- replicate(n,2)
  } else {
    ordinal_levels <- replicate(n,sample(c(2:(2 * exp_levels - 2)),1))
  }
  ordinal_data <- scaled_data
  cuts_list <- vector("list", n)  # List to store cuts for each variable
  for (i in c(1:n)) {
    
    check_levels <- ordinal_levels[i] - 1
    while (check_levels != ordinal_levels[i]) {
      cuts <- c(-Inf,
                cutfun(ordinal_levels[i],concent_param),
                Inf)
      temp <- cut(scaled_data[,i], simplify2array(cuts), labels = FALSE) - 1
      check_levels <- length(unique(temp))
    }
    ordinal_data[,i] <- temp
    cuts_list[[i]] <- cuts
  }
  colnames(ordinal_data) <- c(1:n)
  return(list(ordinal_data = ordinal_data, cuts = cuts_list))
}
  
##' getCov(S, DAG):
##' a function that determines the Cholesky factorization of a covariance matrix in correlation form
##' @param S: covariance matrix in correlation form
##' @param DAG: a DAG estimated using OSEM 
##' @return a list of Cholesky factor B and V   
getCov <- function(S, DAG){
  n <- nrow(DAG)
  B <- matrix(0, nrow = n, ncol = n)

  for (j in c(1:ncol(DAG))) {
    pa <- which(DAG[,j] != 0)
    lp <- length(pa)
    switch (as.character(lp),
          "0" = {#no parents
            next
          },
          "1" = {#one parent
            B[j,pa] <- 1 / S[pa,pa] * S[j,pa]
          },
          {#more parents
            B[j,pa] <- chol2inv(chol(S[pa,pa])) %*% S[j,pa]
          }
    )
  }

  I_B <- diag(n) - B
  V <- I_B%*% S %*% t(I_B)
  
  Chol<-list()
  Chol[[1]]<-B
  Chol[[2]]<-V
  
  return(Chol)
}

# Wrapper function for OwenT to handle infinite inputs
safeOwenT <- function(h, a) {
  if (is.infinite(h)) {
    return(0)  # T(h, a) = 0 if h = ±Inf
  }
  if (is.infinite(a)) {
    # T(h, a) when a = +Inf
    return(0.5 * pnorm(h) * (1 - pnorm(h)))
  }
  if (a == -Inf) {
    return(0)  # T(h, a) = 0 if a = -Inf
  }
  # Call OwenT for finite inputs
  return(OwenT(h, a))
}

##' getT(a,alpha,b):
##' A function that compute pmvnorm(upper=c(alpha, a/sqrt(1+b^2)), sigma=cbind(c(1, -b/sqrt(1+b^2)),c( -b/sqrt(1+b^2),1)))-1/2*pnorm(alpha) 
##' with Owen functions
##' @param a: param of BVN integration
##' @param alpha: param of BVN integration
##' @param b: param in BVN correlation 
##' @return a number representing pmvnorm(upper=c(alpha, a/sqrt(1+b^2)), sigma=cbind(c(1, -b/sqrt(1+b^2)),c( -b/sqrt(1+b^2),1)))-1/2*pnorm(alpha)   
getT <- function(a, alpha, b) {
  rho <- sqrt(1 + b^2)
  
  if (!is.finite(a)) {
    if (alpha == a) {
      return(0)
    } else if (alpha == -a) {
      return(0.5)
    } else {
      return(as.numeric(a / alpha < 0) / 2 + safeOwenT(alpha, a / alpha))
    }
  }
  
  if (!is.finite(alpha)) {
    sgn <- ifelse(a == 0, +1, sign(a))
    return(as.numeric(sgn != sign(alpha)) / 2 + safeOwenT(a / rho, b))
  }
  
  if (a == 0 && alpha == 0) {
    return(safeOwenT(0, 1 + b) + safeOwenT(0, b + rho^2))
  }
  
  if (a == 0) {
    return(safeOwenT(0, alpha + b))
  }
  
  if (alpha == 0) {
    return(safeOwenT(a / rho, b))
  }
  
  sgn <- sign(a)
  sgnal <- sign(alpha)
  return(as.numeric(sgn != sgnal) / 2 + safeOwenT(alpha, a / alpha + b) + safeOwenT(a / rho, b + ((alpha * rho^2) / a)))
}


##' getProb(sigma_i, mean_io, Woi, sigma_do, lev_out, cuts_int, intType):
##' a function that determine post-intervention probability of an outcome level for each intervention level 
##' @param sigma_i: sd of intervention variable
##' @param mean_io: vector (mu_i,mu_o)
##' @param Woi: path effect of intervention variable on outcome variable 
##' @param sigma_do: sd of post-intervention distribution of outcome variable
##' @param out_band: level band of outcome variable
##' @param cuts_int: list of thresholds of intervention variable
##' @param intType: type of integration to be applied ("DIS", "BVN", "OWEN")
##' @return a vector representing probability of lev_out for each intervention level of X_i 
getProb<-function(sigma_i, mean_io, Woi, sigma_do, out_band, cuts_int, intType=c("DIS", "BVN", "OWEN"))
{
  prob<-rep(0,length(cuts_int)-1)
  
  if (intType == "DIS"){
    InnerIntegrald = function(y){pnorm(out_band[2],mean=mean_io[2]+Woi*(y-mean_io[1]), sd=sigma_do)-pnorm(out_band[1],mean=mean_io[2]+Woi*(y-mean_io[1]), sd=sigma_do)}
    ##Alternatively using integrate() 
    #InnerFuncd = function(x, y){dnorm(x, mean = mean_io[2]+Woi*(y-mean_io[1]), sd=sigma_do)}
    #InnerIntegrald = function(y){sapply(y,function(z) {integrate(InnerFuncd, out_band[1],out_band[2], y = z)$value})}
    for (l in seq(2, length(cuts_int))){
      OuterFuncd = function(y){InnerIntegrald(y)*(dtruncnorm(y, mean =mean_io[1], sd = sigma_i, a=cuts_int[l-1],b=cuts_int[l]))}
      prob[l-1]<- round(integrate(OuterFuncd,cuts_int[l-1], cuts_int[l])$value, 6)
    }
  } else {
    #components of standardized outcome band 
    a<- (out_band - mean_io[2])/sigma_do
    b<- -(Woi*sigma_i)/sigma_do
    # standardized intervention band 
    alpha_bar <- (cuts_int - mean_io[1])/sigma_i
    
    if (intType == "BVN"){
      corr<- -b/sqrt(1+b^2)
      #c(BN(alpha_bar(l-1), a[1]), BN(alpha_bar(l-1), a[2]))
      Bnlow<-c(pmvnorm(upper=c(alpha_bar[1], a[1]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE),pmvnorm(upper=c(alpha_bar[1], a[2]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE)) 
      #c(BN(alpha_bar(l), a[1]), BN(alpha_bar(l), a[2]))
      Bnup<- c(pmvnorm(upper=c(alpha_bar[2], a[1]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE),pmvnorm(upper=c(alpha_bar[2], a[2]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE)) 
      philow<-pnorm(alpha_bar[1]) #Phi(alpha_bar(l-1))
      phiup<- pnorm(alpha_bar[2])  #Phi(alpha_bar(l))
      for (l in seq(1, length(cuts_int)-1)){
        prob[l]<- round((Bnup[2]-Bnlow[2]-Bnup[1]+Bnlow[1])/(phiup-philow),6)
        
        if (l == length(cuts_int)-1){
          break 
        }
        Bnlow<-Bnup
        Bnup<- c(pmvnorm(upper=c(alpha_bar[l+2], a[1]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE),pmvnorm(upper=c(alpha_bar[l+2], a[2]/sqrt(1+b^2)), sigma=cbind(c(1,corr),c(corr,1)),keepAttr = FALSE)) 
        philow<- phiup
        phiup<-pnorm(alpha_bar[l+2])
      }
    }
    
    if (intType == "OWEN"){
      #c(T(param(a[1],alpha_bar[l-1]), T(param(a[2], alpha_bar[l-1])))
      Tlow<-c(getT(a[1],alpha_bar[1],b), getT(a[2],alpha_bar[1],b))
      #c(T(param(a[1],alpha_bar[l]), T(param(a[2], alpha_bar[l])))
      Tup<- c(getT(a[1],alpha_bar[2],b), getT(a[2],alpha_bar[2],b))
      philow<-pnorm(alpha_bar[1]) #Phi(alpha_bar(l-1))
      phiup<- pnorm(alpha_bar[2])  #Phi(alpha_bar(l))
      for (l in seq(1, length(cuts_int)-1)){
        prob[l]<- round((-Tup[2]+Tlow[2]+Tup[1]-Tlow[1])/(phiup-philow),6)
        
        if (l == length(cuts_int)-1){
          break 
        }
        Tlow<-Tup
        philow<- phiup
        phiup<-pnorm(alpha_bar[l+2])
        
        Tup<- c(getT(a[1],alpha_bar[l+2],b), getT(a[2],alpha_bar[l+2],b))
      }
    }
  }
  
  return(prob) 
}


##' getEffects(i,mu, B, V, cuts, intType):
##' a function that determine effects of intervention variable X_i on the other variables 
##' @param i: intervention variable 
##' @param mu: mean vector 
##' @param B: path coefficients matrix
##' @param V: noise variance matrix
##' @param cuts: list of thresholds
##' @param intType: type of integration to be applied ("DIS", "BVN", "OWEN")
##' @return a list of matrix representing effects of X_i on each level of the other variables    
getEffects<- function(i,mu, B,V,cuts, intType=c("DIS", "BVN", "OWEN")){
  n <- nrow(B)
  S <- ginv(diag(n)-B) %*% V %*% t(ginv(diag(n)-B))
  sigma_i<- sqrt(S[i,i])
  V_new <- V
  B_new <- B
  V_new[i,i]<-0
  B_new[i,]<-0
  W_new <- solve(diag(nrow(B)) - B_new)
  
  Sigma_do <- W_new%*%V_new%*%t(W_new)
  
  effect_i<-list()
  names<-c()
  nlev_i<- length(cuts[[i]])-1 
  for (o in seq(n)){
    names[o]<-paste0("Ord_",o)
    if (o == i){#outcome variable must be different from intervention variable
      effect_i[[o]]<-array(0,0) #empty 
      next
    }
  Woi<- W_new[o,i] 
  sigma_do <- sqrt(Sigma_do[o,o])
  nlev_o<- length(cuts[[o]])-1
  effect_i[[o]]<-array(0, c(nlev_i,nlev_i,nlev_o))
  for (k in seq(1,nlev_o)){
    temp<-rep(0,nlev_i) #matrix storing on the columns intervention probabilities for each outcome level
    band<-c(cuts[[o]][k], cuts[[o]][k+1])
    mean_io <- c(mu[i], mu[o])
    temp<-getProb(sigma_i, mean_io, Woi, sigma_do, band, cuts[[i]], intType=intType)
    effect_i[[o]][ , ,k]<-t(outer(temp, temp, '-')) #effect on level k of outcome of shifting X_i from level in row to level in col
    # option for rounding
    effect_i[[o]][ , ,k]<-round(effect_i[[o]][ , ,k],4)
  }
}
names(effect_i)<-names
return(effect_i)
}


getallEffects<- function(mu, B,V,cuts, intType=c("DIS", "BVN", "OWEN")){
  n<-nrow(B)
  names<-c()
  effect<-list()
  for (i in c(1:n)){
    names[i]<-paste0("Int_",i)
    effect[[i]]<-getEffects(i,mu, B,V,cuts, intType=intType)
  }
  names(effect)<-names
  return(effect)
}




##' getParam(n, AM, data, usrpar):
##' a function that implements the structural EM algorithm for learning 
##' just the parameters of Bayesian networks from ordinal data given the graph structure (adapted from ordinalStructEM Luo et. all )
##' @param n: dimension of the variable
##' @param AM: Adjacency matrix of the graph structure 
##' @param data: ordinal dataset
##' @param usrpar: list of parameters for ordinal scores
##' @param computeObservedLL: compute the observed likelihood of the parameters (default: FALSE)
##' @param iterMCMC_alpha: significance level for the constraint-based initialization (default: 0.05)
##' @param regsubsets: a boolean value indicating whether subset selection is used in parameter update (default: FALSE)
##' @return an object of class MCMCtrace

getParam <- function(n, AM, data,
                            usrpar = list(penType = c("AIC","BIC","other"),
                                          L = 5,
                                          lambda = 0,
                                          preLevels = NULL),
                            computeObservedLL = FALSE,
                            iterMCMC_alpha = 0.05,
                            regsubsets = FALSE) {
  
  start_time = Sys.time()
  print("Initializing parameters and DAG...")
  param <- scoreparameters("usr",data,usrpar = usrpar)
  
  print("SEM iteration starts...")
  iter <- 0
  currentBestDAG <- AM
  results<- list()
  
  # Maximum iterations to control runtime (can change)
  if ((param$lambda <= 1) || (n >= 30)) {
    nr_EM_iter <- 5
    nr_plus1it <- 5
  } else {
    nr_EM_iter <- 10
    nr_plus1it <- 10
  }
  
  while (iter < nr_EM_iter) {
    
    # Compute expected statistics
    param <- getExpectedStats(param)
    
    # Structure update
    candidateBestDAG <- AM
    
    # Parameter update
    param <- ordinalUpdateParam(param,candidateBestDAG, regsubsets = regsubsets)
    
    currentBestDAG <- candidateBestDAG
    
    iter <- iter + 1
  }
  print("SEM iteration ends...")
  
  if (computeObservedLL) {
    print("Calculating observed log-likelihood...")
    results$observed_score <- observedLL(param) - param$lambda * log(param$N) / 2 * sum(currentBestDAG)
  }
  
  results$param <- param
  
  end_time = Sys.time()
  print(end_time - start_time)
  results$runtime <- as.double(end_time - start_time,units = "secs")
  
  return(results)
  
}



















