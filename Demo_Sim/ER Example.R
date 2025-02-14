########################################################
# Ordinal Causal Effect Estimation on ER DAG for Regenerated Data #
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
library(data.table)


# # Function to count the number of undirected edges in a CPDAG
# count_undirected_edges <- function(cpdag) {
#   adj_matrix <- as(cpdag, "matrix")  # Convert CPDAG to adjacency matrix
#   
#   # Check for symmetric edges (undirected edges)
#   undirected_edges <- sum(adj_matrix != 0 & adj_matrix == t(adj_matrix)) / 2
#   
#   return(undirected_edges)
# }
# 
# # Define function to find a DAG whose CPDAG has no undirected edges
# find_unique_cpdag_dag <- function(n, d, method, wFUN, max_attempts = 10000) {
#   found <- FALSE
#   attempts <- 0
#   seed <- NA
#   target_dag <- NULL
#   
#   while (!found && attempts < max_attempts) {
#     attempts <- attempts + 1
#     set.seed(attempts)  # Use the attempt number as the seed
#     current_seed <- .Random.seed
#     
#     # Generate a random DAG
#     dag <- randDAG(n = n, d = d, method = method, wFUN = wFUN)
#     
#     # Compute the CPDAG
#     cpdag <- dag2cpdag(dag)
#     
#     # Count undirected edges in the CPDAG
#     undirected_edge_count <- count_undirected_edges(cpdag)
#     
#     # If no undirected edges, the equivalence class has only one DAG
#     if (undirected_edge_count == 0) {
#       found <- TRUE
#       seed <- current_seed
#       target_dag <- dag
#     }
#   }
#   
#   if (found) {
#     list(dag = target_dag, seed = seed)
#   } else {
#     stop("No such DAG found within the maximum attempts.")
#   }
# }



### 2. Generate DAGs and ordinal data
n <- 16

# # # NOT RUN
# # # Generate a regular DAG (only element of its equivalence class) with n nodes with 5 number of neighbors
# resultDAG <- find_unique_cpdag_dag(
#   n = n,
#   d = 5,
#   method = "er",
#   wFUN = list(mywFUN)
# )
# 
# trueDAG<- resultDAG$dag
# trueAM <- as(trueDAG,"matrix")
# # Keep only the structure and no weights in AM
# trueAM[which(trueAM != 0)]<-1 
# 
# trueB<- t(as(trueDAG,"matrix"))
# truecov <- trueCov(trueDAG)
# # Transform covariance matrix in its correlation form
# D <- diag(sqrt(diag(truecov)))
# D.inv <- chol2inv(chol(D))
# trueSigma <- D.inv %*% truecov %*% D.inv
# 
# save(trueDAG,n,trueAM, trueB, truecov, D, D.inv,trueSigma, resultDAG, file = "./Demo_Sim/ERexample_data.RData")
load("Demo_Sim/ERexample_data.RData")

# Plot True Dag 
pdf("TrueDAG.pdf", width = 11.55, height = 10.69)
plotDAG<- plot(trueDAG)
dev.off()

### 3. Generate Gaussian dataset and convert into an ordinal dataset
# Size of the Dataset 
N <- 500 

# # Get the Latent Datasets 
# hidden_data <- rmvDAG2(N, trueDAG)
# scaled_data <- t(t(hidden_data) - apply(hidden_data,2,mean))
# True_scaled_data <- t(D.inv %*% t(scaled_data))
# 
# # Extract cuts and generate True Ordinal 
# Trueres<-extractCuts(True_scaled_data, exp_levels = 4, concent_param = 2)
# TrueCuts <- Trueres$cuts
# TrueOrdinal <-Trueres$ordinal_data
# 
# ### 4. Compute True Ordinal Causal Effects 
# CholTrue <- getCov(trueSigma,trueAM)
# Beff <- CholTrue[[1]]
# Veff <- CholTrue[[2]]
# TrueEffects<- getallEffects(mu=rep(0,n),B=Beff,V=Veff,TrueCuts,intType = "DIS")
# TrueEffects_BVN<- getallEffects(mu=rep(0,n),B=Beff,V=Veff,TrueCuts,intType = "BVN")
# # # They are all equal Sanity Check
# identical(TrueEffects, TrueEffects_BVN)
# 
# ### 5. Generate Simulation Sample
# 
# boot.sample<-list()
# for (i in 1:500){
#   hidden_data <- rmvDAG2(N, trueDAG)
#   scaled_data <- t(t(hidden_data) - apply(hidden_data,2,mean))
#   scaled_data <- t(D.inv %*% t(scaled_data))
#   boot.sample[[i]]<- convertcutsToOrdinal(scaled_data, cuts=TrueCuts)
# }
# 
# save(True_scaled_data, TrueCuts, TrueOrdinal, boot.sample, file = "./Demo_Sim/BootER_List.RData")
# save(CholTrue, Beff, Veff, TrueEffects, file = "./Demo_Sim/TrueEffects.RData")
load("Demo_Sim/BootER_List.RData")
load("Demo_Sim/TrueEffects.RData")



### 6. Learn Parameters
# NOT RUN
# source("ER_example.boot.R")
load("./Demo_Sim/Bootresults_ER.RData")

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
# 
# 
# 
# #### 8. Generate list of effects for each int-out couple 
# 
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
# save(effects4couple_OSEM, file = "./Demo_Sim/ER_Effects_OSEM.RData")
# save(effects4couple_param, file = "./Demo_Sim/ER_Effects_param.RData")
# load("Demo_Sim/ER_Effects_OSEM.RData")
# load("Demo_Sim/ER_Effects_param.RData")



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
# Demo_data <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
# 
# # Convert to appropriate types
# Demo_data[, `:=` (
#   Method = factor(Method, levels = c("BN", "Param")),
#   Int = factor(Int, levels = sort(unique(Int)), ordered = TRUE),
#   MaxLevelInt = as.integer(MaxLevelInt),
#   Outcome = factor(Outcome, levels = sort(unique(Outcome)), ordered = TRUE),
#   Level = factor(Level, levels = sort(unique(Level)), ordered = TRUE),
#   OCE = as.numeric(OCE)
# )]

# save(Demo_data, file ="Demo_Sim/Demo_Data.RData")
load("Demo_Sim/Demo_Data.RData")


















# Compute variance of OCE grouped by Int, Outcome, Level, and Method
Demo_oce_variance <- Demo_data[, .(Variance_OCE = var(OCE, na.rm = TRUE)), by = .(Int, Outcome, Level, Method)]

# Print the result
print(Demo_oce_variance)


# Comparison of variances between the approaches  
# Rename variance columns before merging to avoid conflicts
setnames(Demo_oce_variance, "Variance_OCE", "Variance_OCE_Demo")
setnames(Boot_oce_variance, "Variance_OCE", "Variance_OCE_Boot")
# Merge datasets on common keys (Int, Outcome, Level, Method)
oce_comparison <- merge(
  Demo_oce_variance, Boot_oce_variance, 
  by = c("Int", "Outcome", "Level", "Method"), 
  suffixes = c("Demo", "Boot")
)


# Compute absolute and relative differences in variance
oce_comparison[, `:=` (
  Variance_Diff = Variance_OCE_Demo - Variance_OCE_Boot,  # Absolute Difference
  Variance_Ratio = Variance_OCE_Demo / Variance_OCE_Boot  # Ratio (Relative Difference)
)]

# Print results
print(oce_comparison)



# # Compute MSE Error
# 
Sigma_Param<- list()
Sigma_OSEM<- list()
Blist<-list()
Vlist <-list()
BlistO<-list()
VlistO <-list()
cutslist<-list()
cuts0<-list()

for (i in 1:500){
  S<- results[[i]]$Param$S
  Sigma_Param[[i]]<-S
  Chol <- getCov(S,trueAM)
  Blist[[i]] <- Chol[[1]]
  Vlist[[i]] <- Chol[[2]]
  cutslist[[i]]<-results[[i]]$Param$cuts

  S<- results[[i]]$OSEM$S
  DAG <- as.matrix(results[[i]]$OSEM$DAG)
  Sigma_OSEM[[i]]<-S
  Chol <- getCov(S,DAG)
  BlistO[[i]] <- Chol[[1]]
  VlistO[[i]] <- Chol[[2]]
  cuts0[[i]]<-results[[i]]$OSEM$cuts
}

# Function to compute MSE for a matrix
compute_mse_matrix <- function(estimated_matrix, true_matrix) {
  return(mean((estimated_matrix - true_matrix)^2))
}



# Initialize lists to store element-wise MSE matrices
mse_matrixParam_list <- list()
mse_BPar_list <- list()
mse_matrixOSEM_list <- list()
mse_BO_list <- list()

# Function to compute element-wise MSE for matrices
compute_elementwise_mse <- function(estimated_matrix, true_matrix) {
  return((estimated_matrix - true_matrix)^2)  # Element-wise squared error
}

# Compute element-wise MSE for each simulation
for (i in 1:500) {
  mse_matrixParam_list[[i]] <- compute_elementwise_mse(Sigma_Param[[i]], trueSigma)
  mse_BPar_list[[i]] <- compute_elementwise_mse(Blist[[i]], Beff)
  
  mse_matrixOSEM_list[[i]] <- compute_elementwise_mse(Sigma_OSEM[[i]], trueSigma)
  mse_BO_list[[i]] <- compute_elementwise_mse(BlistO[[i]], Beff)
}

# Convert lists to arrays for easier aggregation (assuming all matrices have the same dimensions)
mse_matrixParam_array <- simplify2array(mse_matrixParam_list)  # 3D array (rows, cols, simulations)
mse_BPar_array <- simplify2array(mse_BPar_list)
mse_matrixOSEM_array <- simplify2array(mse_matrixOSEM_list)
mse_BO_array <- simplify2array(mse_BO_list)

# Compute mean MSE per element across simulations
mean_mse_matrixParam <- apply(mse_matrixParam_array, c(1,2), mean, na.rm = TRUE)
mean_mse_BPar <- apply(mse_BPar_array, c(1,2), mean, na.rm = TRUE)
mean_mse_matrixOSEM <- apply(mse_matrixOSEM_array, c(1,2), mean, na.rm = TRUE)
mean_mse_BO <- apply(mse_BO_array, c(1,2), mean, na.rm = TRUE)

# Print or inspect results
print(mean_mse_matrixParam)
print(mean_mse_BPar)
print(mean_mse_matrixOSEM)
print(mean_mse_BO)

# Aggregate MSE across the simulations 
# Compute MSE for each simulation
mse_matrixParam <- sapply(1:500, function(i) compute_mse_matrix(Sigma_Param[[i]], trueSigma))
mse_BPar <- sapply(1:500, function(i) compute_mse_matrix(Blist[[i]], Beff))

mse_matrixOSEM <- sapply(1:500, function(i) compute_mse_matrix(Sigma_OSEM[[i]], trueSigma))
mse_BO <- sapply(1:500, function(i) compute_mse_matrix(BlistO[[i]], Beff))



matrix_dim <- dim(trueSigma)  # Get the dimension of the matrix
# Initialize an array to store MSE values for each matrix
mse_array <- array(0, dim = c(matrix_dim[1], matrix_dim[2], 500))


# Compute element-wise MSE for each matrix
for (i in 1:500) {
  mse_array[,,i] <- compute_elementwise_mse(Sigma_OSEM[[i]], trueSigma)
}

# Compute the mean across all 500 matrices for each element
elementwise_mse <- apply(mse_array, c(1, 2), mean)

# Print or visualize the element-wise MSE
print(elementwise_mse)

# Optional: Heatmap visualization (if you want to see a heatmap of MSE)
heatmap(elementwise_mse, main = "Element-wise MSE", col = heat.colors(100))

library(gplots)

# Generate heatmap with legend
heatmap.2(elementwise_mse,
                      main = "Element-wise MSE True vs OSEM Sigma",
                        trace = "none",
                        col = heat.colors(100),
                         key = TRUE,   # Enables the color legend
                         key.title = "MSE",  # Title for the color legend
                         )  # Hides density plot in legend

