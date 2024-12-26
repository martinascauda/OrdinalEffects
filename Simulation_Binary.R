library(mnormt)
library(MASS)
library(truncnorm)
# Mean vector
mu <- rep(0,2)
mu_1 = mu[1]
mu_2 = mu[2]

# Covariance matrix
b <- 0.5
B <- matrix(c(0,0, b,0), nrow=2,ncol=2,byrow= TRUE)
v<-rep(1,2)
v_1=v[1]
v_2 = v[2]
V <- matrix(c(v_1,0, 0,v_2), nrow=2,ncol=2,byrow= TRUE)
I <- diag(2)
# Covariance matrix
S <- ginv(I-B) %*% V %*% t(ginv(I-B))
s_1 <- S[1,1]
s_2 <- S[2,2]

# Thresholds
alpha<-c(0.2, 0.4)
alpha_1 <- alpha[1]
alpha_2 <- alpha[2]

# Contour level plot
xsim <- seq(-3, 5, 0.1)
ysim<-seq(-5, 3, 0.1)
f <- function(xsim,ysim) dmnorm(cbind(xsim,ysim), mu, S)
z <- outer(xsim,ysim,f)
filled.contour(xsim,ysim,z, plot.title={
  title(main = "Density estimation: contour plot", xlab=expression("Y"[1]), ylab=expression("Y"[2]))
  abline(v=alpha[1], col = "red", lwd = 1)
  abline(h=alpha[2], col = "red", lwd = 1)
})

# BINARY CASE - CLASSIC APPROACH
set.seed(123)

# Probabilities(sum=1)
tau<- sadmvn(lower=rep(-Inf,2), upper=alpha, mean=mu, varcov =S)
beta<- sadmvn(lower=c(-Inf,alpha[2]), upper=c(alpha[1],Inf), mean=mu, varcov =S)
gamma<- sadmvn(lower=c(alpha[1],-Inf), upper=c(Inf,alpha[2]), mean=mu, varcov =S)
theta<- sadmvn(lower=alpha, upper=rep(Inf,2), mean=mu, varcov =S)

# Intervention effects
I0 = (gamma)/(gamma+theta) - (tau)/(tau+beta)
I1 = (theta)/(theta+gamma)-(beta)/(tau+beta)

# BINARY CASE - LATENT APPROACH
# x is y2 while y is y1 in the toy model

# Strategy 1 - Distribution approach 
InnerFuncd = function(x, y){dnorm(x, mean = mu_2+b*(y-mu_1), sd=sqrt(v_2))*(dtruncnorm(y, mean =mu_1, sd = sqrt(s_1), a=alpha_1,b=+Inf))}
InnerFuncdprime = function(x, y){dnorm(x, mean = mu_2+b*(y-mu_1), sd=sqrt(v_2))*(dtruncnorm(y, mean =mu_1, sd = sqrt(s_1), a=-Inf,b=alpha_1))}
InnerIntegrald = function(x){sapply(x,function(z) {integrate(InnerFuncd, alpha_1, Inf, x = z)$value})}
InnerIntegraldprime = function(x){sapply(x,function(z) {integrate(InnerFuncdprime, -Inf,alpha_1, x = z)$value})}
I0d <- integrate(function(x){InnerIntegrald(x) - InnerIntegraldprime(x)}, -Inf, alpha_2)$value
I1d<- integrate(function(x){InnerIntegrald(x) - InnerIntegraldprime(x)}, alpha_2, Inf)$value

# Strategy 2 - Quantile approach
InnerFuncq = function(x, y){(dnorm(x, mean = mu_2+b*(y-mu_1), sd=sqrt(v_2))-
                              dnorm(x, mean = mu_2+b*(qtruncnorm(ptruncnorm(y, mean =mu_1, sd = sqrt(s_1), a=alpha_1,b=+Inf), a=-Inf, b=alpha_1, mean =mu_1, sd = sqrt(s_1))-mu_1), 
                                    sd=sqrt(v_2)))*dtruncnorm(y, mean =mu_1, sd = sqrt(s_1), a=alpha_1,b=+Inf)}
InnerIntegralq = function(x){sapply(x,function(z) {integrate(InnerFuncq, alpha_1, Inf, x = z)$value})}
I0q <- integrate(InnerIntegralq, -Inf, alpha_2)$value
I1q<- integrate(InnerIntegralq, alpha_2, Inf)$value

# Print results 
paste("I0 =", I0, "I0distribution=", I0d, "I0quantile=", I0q)
paste("I1 =", I1, "I1distribution=", I1d, "I1quantile=", I1q)




