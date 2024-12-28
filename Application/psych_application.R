########### Psychological Application ###########
# Required packages
library(BiDAG)
library(pcalg)
library(bnlearn)
library(corrplot)
library(igraph)
library(abind)
source("ordinalEffects.R")
source("./OSEMSource/R/OrdinalScore.R")
insertSource("./OSEMSource/R/spacefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/usrscorefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/initpar.R",package = "BiDAG")
insertSource("./OSEMSource/R/scoreagainstdag.R",package = "BiDAG")

# Data processing
load("./Application/OCDRogers.RData")
n <- ncol(datRogers)
N <- nrow(datRogers)

load("./Application/Results_boot_pysch.RData")

###########
# Bootstrapping
# NOT RUN
#source("psych_boot.R")

# Compute Effects
# NOT RUN

#OWEN Methods
# Effects<-list()
# system.time(for (i in 1:500){
# OSEMfit<-results[[i]]
#   cuts <-OSEMfit[["param"]][["cuts"]]
#   S<- OSEMfit[["param"]][["Sigma_hat"]]
#   DAG <- as.matrix(OSEMfit$DAG)
#   Chol <- getCov(S,DAG)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effects[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "OWEN")
#})

## Distributional method
# Effectsdis<-list()
# system.time(for (i in 1:500){
#   OSEMfit<-results[[i]]
#   cuts <-OSEMfit[["param"]][["cuts"]]
#   S<- OSEMfit[["param"]][["Sigma_hat"]]
#   DAG <- as.matrix(OSEMfit$DAG)
#   Chol <- getCov(S,DAG)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effectsdis[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "DIS")
#   print(i)
# })

#save Effects to a file 
#save(Effects, file = "Psych_Effects.RData")
# save(Effectsdis, file="Psych_EffectsDIS.RData")
load("./Application/Psych_Effects.RData")
#load("./Application/Psych_EffectsDIS.RData")

# Check computations
# for (i in 1:500){
#   if (!(identical(round(unlist(Effects[[i]]),2),round(unlist(Effectsdis[[i]]),2)))){
#   print(i)
#   }
# }

# i = 338 not equal, check the position of the mismatch (code below): int 5 on out 8. 
# for (k in 1:24){
#   for (j in 1:24){
#   if (!(identical(round(unlist(Effects[[i]][[k]][[j]]),2),round(unlist(Effectsdis[[i]][[k]][[j]]),2)))){
#     print(paste0(k,j))
#   }
#   }
# }



#Generate list of effects for each int-out couple
effects4couple<-vector(mode="list", 24)
for (i in 1:24){
  effects4couple[[i]]<-vector(mode="list",24)
  for (j in 1:24){
    effects4couple[[i]][[j]]<-vector(mode="list",500)
    for (k in c(1:500)){
      effects4couple[[i]][[j]][[k]]<-Effects[[k]][[i]][[j]]
}
}
}


# Point estimates
OSEMfit_point <- ordinalStructEM(n, datRogers,
                           usrpar = list(penType = "other",
                                         L = 5,
                                         lambda = 6))

g <- as_graphnel(graph_from_adjacency_matrix(OSEMfit_point$DAG))
cpdag_OSEM <- dag2cpdag(g)
png(paste0("OSEM_CDDAG_Pysch",i,".png"), width = 465, height = 225, units='mm', res = 300)
plot(as(cpdag_OSEM, "graphNEL"),main = "Cpdag estimated with OSEM")
dev.off()


# Get strengths arrows 
#Credits: OSEM code source application 
# res<-list()
# for (i in 1:500){
# res[[i]]<- results[[i]]$maxtrace[[1]]$DAG
# }
# newarray<-array(NA, dim=c(n,n,500))
# for (i in 1:n){
#   for (j in 1:n){
#     for (k in 1:500){
#     newarray[i,j,k]<-res[[k]][i,j]
#     }
#   }
# }
# 
# res_cpdag<-newarray
# for (j in c(1:500)) {
#   res_cpdag[,,j] <- dag2cpdag(newarray[,,j])
# }
# res_OSEM<-res_cpdag
# save(res_OSEM, file="Res_OSEM_Pysch.RData")
load("./Application/Res_OSEM_Pysch.RData")

get_strength <- function (c, cboot) {
  n <- nrow(c)
  for (i in c(1:n)) {
    for (j in c(1:n)) {
      if (c[i,j] & i != j) {
        c[i,j] <- mean(apply(cboot,3,function (A) if ((A[i,j] + A[j,i]) == 1) {A[i,j]} else {A[i,j] / 2}))
      }
    }
  }
  return(c)
}
cpdag_OSEM<- as(OSEMfit_point$DAG,"matrix")
cpdag_OSEM <- get_strength(cpdag_OSEM, res_OSEM)
corrplot(cpdag_OSEM, method = "shade", is.corr = FALSE,
         tl.col = "grey", col.lim = c(0,1),
         mar = c(1,0,0,0)+0.5,addgrid.col = "lightgrey", diag=FALSE)
title(sub = "OSEM")





######### Plots Rain PLot New 
# Preallocate the data frame for plotting Effects 9 on 10 
Effects910 <- effects4couple[[9]][[10]]
total_rows <- 500 * 4 * 6  # Total iterations
data <- data.frame(Level = integer(total_rows),
                   Change = integer(total_rows),
                   OCE = numeric(total_rows))
 
row_index <- 1
for (i in 1:500) {
  for (l in 1:4) {
    for (k in 1:6) {
      if (k %in% c(1, 2, 3)) {
        data[row_index, ] <- c(l, k, Effects910[[i]][,,l][1, k + 1])
      } else if (k %in% c(4, 5)) {
        data[row_index, ] <- c(l, k, Effects910[[i]][,,l][2, k - 1])
      } else {
        data[row_index, ] <- c(l, k, Effects910[[i]][,,l][3, 4])
      }
      row_index <- row_index + 1
    }
  }
}
data$Change <- as.factor(data$Change)
data$Level <- as.factor(data$Level)

save(data, file = "Application/DataInt910_For_Plots.RData")







