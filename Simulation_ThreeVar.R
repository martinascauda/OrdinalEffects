library(mnormt)
library(MASS)
library(truncnorm)

# Mean vector
mu <- rep(0,3)
mu_1 = mu[1]
mu_2 = mu[2]
mu_3 = mu[3]

# Covariance matrix
b12 <- 0.5
b13 <- 0.8
b23 <- 0.9
B <- matrix(c(0,0,0,b12,0,0,b13,b23,0), nrow=3,ncol=3,byrow= TRUE)

v_1 = 1
v_2 = 1
v_3 = 1
V <- matrix(c(v_1,0,0, 0,v_2,0,0,0,v_3), nrow=3,ncol=3,byrow= TRUE)
I <- diag(3)
# Covariance matrix
S <- ginv(I-B) %*% V %*% t(ginv(I-B))
s_1 <- S[1,1]
s_2 <- S[2,2]
s_3 <- S[3,3]

# Thresholds
alpha<-c(1.2, 2.4, 3.3)
alpha_1 <- alpha[1]
alpha_2 <- alpha[2]
alpha_3 <- alpha[3]

# BINARY CASE - CLASSIC APPROACH
set.seed(123)

# Probabilities(sum=1) #remember to take away attributes
p1<- sadmvn(lower=c(alpha_1,-Inf,-Inf), upper=rep(Inf,3), mean=mu, varcov =S) #P(X1 = 1)
p2<- sadmvn(lower=c(-Inf,alpha_2,-Inf), upper=c(alpha_1,Inf,Inf), mean=mu, varcov =S)/(1-p1)#P(X2=1|X1=0)
p3 <- sadmvn(lower=c(alpha_1,alpha_2,-Inf), upper=rep(Inf,3), mean=mu, varcov =S)/(p1)#P(X2=1|X1=1)
tau<-sadmvn(lower=rep(-Inf,3), upper=c(alpha_1,alpha_2,Inf), mean=mu, varcov =S) #P(X1=0,X2=0)
gamma <- sadmvn(lower=c(alpha_1,-Inf, -Inf), upper=c(Inf,alpha_2,Inf), mean=mu, varcov =S) #P(X1=1,X2=0)
beta <- sadmvn(lower=c(-Inf,alpha_2, -Inf), upper=c(alpha_1,Inf,Inf), mean=mu, varcov =S) #P(X1=0,X2=1)
theta <- sadmvn(lower=c(alpha_1,alpha_2, -Inf), upper=rep(Inf,3), mean=mu, varcov =S) #P(X1=1,X2=1)
p4<- sadmvn(lower=c(-Inf,-Inf,alpha_3), upper=c(alpha_1,alpha_2,Inf), mean=mu, varcov =S)/(tau)#P(X3=1|X1=0, X2=0)
p5<- sadmvn(lower=c(alpha_1,-Inf,alpha_3), upper=c(Inf,alpha_2,Inf), mean=mu, varcov =S)/(gamma)#P(X3=1|X1=1, X2=0)
p6<- sadmvn(lower=c(-Inf,alpha_2,alpha_3), upper=c(alpha_1,Inf,Inf), mean=mu, varcov =S)/(beta)#P(X3=1|X1=0, X2=1)
p7<-sadmvn(lower=alpha, upper=rep(Inf,3), mean=mu, varcov =S)/(theta)#P(X3=1|X1=1, X2=1)

# Intervention effects

I0real =(1-p6)*(1-p1)+(1-p7)*p1-(1-p4)*(1-p1)-(1-p5)*p1
I1real = (p6)*(1-p1)+(p7)*p1-(p4)*(1-p1)-(p5)*p1

# BINARY CASE - LATENT APPROACH
# x is effect variable while y is intervention variable in the toy model

# Strategy 1 - Distribution approach
I0d1= integrate(Vectorize(function(x) { 
    integrate(Vectorize(function(y) { 
      integrate(function(p) { 
        dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf))## function you want to calculate
      }, -Inf, Inf)$value     ## 2nd level interval
    }),alpha_2, Inf)$value        ## 1st level interval
  }), -Inf, alpha_3)$value

I0d2= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    integrate(function(p) { 
      dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=-Inf,b=alpha_2))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),- Inf, alpha_2)$value        ## 1st level interval
}), -Inf, alpha_3)$value

I0d = I0d1-I0d2 

I1d1= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    integrate(function(p) { 
      dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),alpha_2, Inf)$value        ## 1st level interval
}), alpha_3, Inf)$value

I1d2= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    integrate(function(p) { 
      dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=-Inf,b=alpha_2))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),- Inf, alpha_2)$value        ## 1st level interval
}), alpha_3, Inf)$value

I1d = I1d1-I1d2 

# Alternative method to compute integral
I0trial1= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    (dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf))*integrate(function(p) { 
      dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),alpha_2, Inf)$value        ## 1st level interval
}), -Inf, alpha_3)$value

I0trial2= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    (dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=-Inf,b=alpha_2))*integrate(function(p) { 
      dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))*dnorm(p, mean = mu_1, sd=sqrt(s_1))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),- Inf, alpha_2)$value        ## 1st level interval
}), -Inf, alpha_3)$value

I0trial=I0trial1-I0trial2

# Strategy 2 - Quantile Approach 

I0q= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    integrate(function(p) { 
      (dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))-dnorm(x, mean = mu_3+b23*(qtruncnorm(ptruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf), a=-Inf, b=alpha_2, mean =mu_2, sd = sqrt(s_2))-mu_2)+b13*(p-mu_1), sd=sqrt(v_2)))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),alpha_2, Inf)$value        ## 1st level interval
}), -Inf, alpha_3)$value

I1q= integrate(Vectorize(function(x) { 
  integrate(Vectorize(function(y) { 
    integrate(function(p) { 
      (dnorm(x, mean = mu_3+b23*(y-mu_2)+b13*(p-mu_1), sd=sqrt(v_2))-dnorm(x, mean = mu_3+b23*(qtruncnorm(ptruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf), a=-Inf, b=alpha_2, mean =mu_2, sd = sqrt(s_2))-mu_2)+b13*(p-mu_1), sd=sqrt(v_2)))*dnorm(p, mean = mu_1, sd=sqrt(s_1))*(dtruncnorm(y, mean =mu_2, sd = sqrt(s_2), a=alpha_2,b=+Inf))## function you want to calculate
    }, -Inf, Inf)$value     ## 2nd level interval
  }),alpha_2, Inf)$value        ## 1st level interval
}), alpha_3, Inf)$value

# Print results 
paste("I0 =", I0real, "I0distribution=", I0d, "I0quantile=", I0q)
paste("I1 =", I1real, "I1distribution=", I1d, "I1quantile=", I1q)







