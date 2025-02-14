########################################################
# Ordinal Causal Effect Estimation on ER DAG for Bootstrapped Data#
######################################################## 


### 1. Load necessary packages and source code
library(Rgraphviz)
library(graph)
library(igraph)
library(pcalg)
source("OrdinalEffects.R")
# Modify some of the existing functions in the BiDAG package to accommodate OSEM user-defined functions
source("./OSEMSource/R/OrdinalScore.R")
insertSource("./OSEMSource/R/spacefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/usrscorefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/initpar.R",package = "BiDAG")
insertSource("./OSEMSource/R/scoreagainstdag.R",package = "BiDAG")
library(BiDAG)
library(dplyr)
library(data.table)

### 2. Generate DAGs and ordinal data
n <- 16
load("Demo_Sim/ERexample_data.RData")

# Plot True Dag 
plotDAG<- plot(trueDAG)

### 3. Generate Gaussian dataset and convert into an ordinal dataset
# Size of the Dataset 
N <- 500 

load("Demo_Sim/BootER_List.RData")
# # ### 5. Generate Simulation Sample
# bootstrap.sample<-list()
# for (i in 1:500){
# bootstrap.sample[[i]]<- TrueOrdinal[sample(N,N,replace = TRUE),]
# }
# save(bootstrap.sample, file = "Boot_Sim/BootstrapER_List.RData")




### 6. Learn Parameters
# NOT RUN
# source("ER_example.boot.R")
# load("Boot_Sim/Bootstrap_results_ER.RData")

### 7. Compute Ordinal Causal Effects 
# NOT RUN
# Effects_OSEM<-list()
# system.time(for (i in 1:500){
#   OSEMfit<-results[[i]]$OSEM
#   cuts <-OSEMfit$cuts
#   S<- OSEMfit$S
#   DAG <- as.matrix(OSEMfit$DAG)
#   Chol <- getCov(S,DAG)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effects_OSEM[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "DIS")
# })
# 
# Effects_param<-list()
# system.time(for (i in 1:500){
#   Paramfit<-results[[i]]$Param
#   cuts <-Paramfit$cuts
#   S<- Paramfit$S
#   Chol <- getCov(S,trueAM)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effects_param[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "DIS")
# })
# # 
# # 
# # 
# # #### 8. Generate list of effects for each int-out couple 
# # 
# effects4couple_OSEM<-vector(mode="list", n)
# for (i in 1:n){
#   effects4couple_OSEM[[i]]<-vector(mode="list",n)
#   for (j in 1:n){
#     effects4couple_OSEM[[i]][[j]]<-vector(mode="list",500)
#     for (k in c(1:500)){
#       effects4couple_OSEM[[i]][[j]][[k]]<-Effects_OSEM[[k]][[i]][[j]]
#     }
#   }
# }
# 
# effects4couple_param<-vector(mode="list", n)
# for (i in 1:n){
#   effects4couple_param[[i]]<-vector(mode="list",n)
#   for (j in 1:n){
#     effects4couple_param[[i]][[j]]<-vector(mode="list",500)
#     for (k in c(1:500)){
#       effects4couple_param[[i]][[j]][[k]]<-Effects_param[[k]][[i]][[j]]
#     }
#   }
# }
# 
# #save Effects to a file
# save(effects4couple_OSEM, file = "./Boot_Sim/ER_Effects_OSEM.RData")
# save(effects4couple_param, file = "./Boot_Sim/ER_Effects_param.RData")
# load("Boot_Sim/ER_Effects_OSEM.RData")
# load("Boot_Sim/ER_Effects_param.RData")


# # Preallocate a list for efficiency
# results_list <- vector("list", n * 500 * (n - 1) * 2)  # Adjust size estimate if necessary
# index <- 1
# 
# # Loop efficiently
# for (int in 1:n) {
#   for (i in 1:500) {
#     for (o in setdiff(1:n, int)) {
# 
#       # Get the max level
#       maxint <- nrow(effects4couple_param[[int]][[o]][[i]][,,1])
#       outlevels <- length(effects4couple_param[[int]][[o]][[i]][1,1,])
# 
#       # Process "Param" and "BN" data together
#       for (l in 1:outlevels) {
#         results_list[[index]] <- list(
#           Method = "Param",
#           Int = int,
#           MaxLevelInt = maxint,
#           Outcome = o,
#           Level = l,
#           OCE = effects4couple_param[[int]][[o]][[i]][,,l][1, maxint]
#         )
#         index <- index + 1
# 
#         results_list[[index]] <- list(
#           Method = "BN",
#           Int = int,
#           MaxLevelInt = maxint,
#           Outcome = o,
#           Level = l,
#           OCE = effects4couple_OSEM[[int]][[o]][[i]][,,l][1, maxint]
#         )
#         index <- index + 1
#       }
#     }
#   }
# }

# # Convert list to data.table in one step (fast)
# Boot_data <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
# 
# # Convert to appropriate types
# Boot_data[, `:=` (
#   Method = factor(Method, levels = c("BN", "Param")),
#   Int = factor(Int, levels = sort(unique(Int)), ordered = TRUE),
#   MaxLevelInt = as.integer(MaxLevelInt),
#   Outcome = factor(Outcome, levels = sort(unique(Outcome)), ordered = TRUE),
#   Level = factor(Level, levels = sort(unique(Level)), ordered = TRUE),
#   OCE = as.numeric(OCE)
# )]
# 
# save(Boot_data, file ="Boot_Sim/Demo_Data.RData")
load("Boot_Sim/Demo_Data.RData")


# Compute variance of OCE grouped by Int, Outcome, Level, and Method
Boot_oce_variance <- Boot_data[, .(Variance_OCE = var(OCE, na.rm = TRUE)), by = .(Int, Outcome, Level, Method)]

# Print the result
print(Boot_oce_variance)









# # Compute MSE Error
# 
# Sigma_Param<- list()
# Sigma_OSEM<- list()
# Blist<-list()
# Vlist <-list()
# BlistO<-list()
# VlistO <-list()
# 
# for (i in 1:500){
#   S<- results[[i]]$Param$S
#   Sigma_Param[[i]]<-S
#   Chol <- getCov(S,trueAM)
#   Blist[[i]] <- Chol[[1]]
#   Vlist[[i]] <- Chol[[2]]
#   
#   S<- results[[i]]$OSEM$S
#   DAG <- as.matrix(results[[i]]$OSEM$DAG)
#   Sigma_OSEM[[i]]<-S
#   Chol <- getCov(S,DAG)
#   BlistO[[i]] <- Chol[[1]]
#   VlistO[[i]] <- Chol[[2]]
# }
# 
# # Function to compute MSE for a matrix
# compute_mse_matrix <- function(estimated_matrix, true_matrix) {
#   return(mean((estimated_matrix - true_matrix)^2))
# }
# 
# # Compute MSE for each simulation
# mse_matrixParam <- sapply(1:500, function(i) compute_mse_matrix(Sigma_Param[[i]], trueSigma))
# mse_BPar <- sapply(1:500, function(i) compute_mse_matrix(Blist[[i]], Beff))
# 
# mse_matrixOSEM <- sapply(1:500, function(i) compute_mse_matrix(Sigma_OSEM[[i]], trueSigma))
# mse_BO <- sapply(1:500, function(i) compute_mse_matrix(BlistO[[i]], Beff))

