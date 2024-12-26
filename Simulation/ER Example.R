########################################################
# Ordinal Causal Effect Estimation on ER DAG with Random cuts          #
######################################################## 


### 1. Load necessary packages and source code
library(Rgraphviz)
library(graph)
library(igraph)
library(pcalg)
source("OrdinalEffects.R")
source("PlotEffects_Rain.R")
# Modify some of the existing functions in the BiDAG package to accommodate our user-defined functions
source("./OSEMSource/R/OrdinalScore.R")
insertSource("./OSEMSource/R/spacefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/usrscorefns.R",package = "BiDAG")
insertSource("./OSEMSource/R/initpar.R",package = "BiDAG")
insertSource("./OSEMSource/R/scoreagainstdag.R",package = "BiDAG")
library(BiDAG)

load("./Simulation/ER_exampleRand.data.RData")


### 2. Generate DAGs and ordinal data
set.seed(123)
N <- 500
n <- 20

# # NOT RUN
# # Generate a regular DAG with 20 nodes with 5 number of neighbors
# trueDAG <- randDAG(n = n, d = 5, method = "er", wFUN = list(mywFUN))
# 
# AM <- as(trueDAG,"matrix")
# AM[which(AM != 0)]<-1 
# 
# B <- t(as(trueDAG,"matrix"))
# mu<-rep(0,nrow(B))
# #variance is unidentifiable, set to 1
# V<-diag(nrow(B))
# 
# save(trueDAG, V,n,N,AM, file = "ER_exampleRand.data.RData")


### 3. Generate Gaussian dataset and convert into an ordinal datasets for the bootstrap
### NOT RUN
# set.seed(12345)
# boot.sample<-list()
# for (i in 1:500){
# boot.sample[[i]]<- generateOrdinal(N,n,trueDAG,exp_levels = 4, concent_param = 2)
# }
# save(boot.sample, file = "BootER_Rand_List.RData")


### 4. Learn parameters using OSEM 

###########
# Bootstrapping
# NOT RUN
# source("ER_example.boot.R")
load("./Simulation/Bootresults_ER_Rand.RData")

# Compute Effects
# NOT RUN

# OSEM (Learning of structure and parameters)
#OWEN Methods
# Effects_OSEM<-list()
# system.time(for (i in 1:500){
#   OSEMfit<-results[[i]]$OSEM
#   cuts <-OSEMfit$cuts
#   S<- OSEMfit$S
#   DAG <- as.matrix(OSEMfit$DAG)
#   Chol <- getCov(S,DAG)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effects_OSEM[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "OWEN")
# })

#Distributional Methods
# Effects_OSEM_dis<-list()
# system.time(for (i in 1:500){
#   OSEMfit<-results[[i]]$OSEM
#   cuts <-OSEMfit$cuts
#   S<- OSEMfit$S
#   DAG <- as.matrix(OSEMfit$DAG)
#   Chol <- getCov(S,DAG)
#   B <- Chol[[1]]
#   Vchol <- Chol[[2]]
#   Effects_OSEM_dis[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "DIS")
# })

#Learning of just parameters
# #OWEN Methods
# Effects_param<-list()
# system.time(for (i in 1:500){
#   Paramfit<-results[[i]]$Param
#   cuts <-Paramfit$cuts
#   S<- Paramfit$S
#   B <- Paramfit$B
#   Vchol <- Paramfit$Vchol
#   Effects_param[[i]] <-getallEffects(mu=rep(0,n),B,V=Vchol,cuts,intType = "OWEN")
# })




#save Effects to a file 
#save(Effects_OSEM, file = "ER_Effects_OSEM.RData")
#save(Effects_OSEM_dis, file = "ER_Effects_OSEM_dis.RData")
#save(Effects_param, file = "ER_Effects_param.RData")
load("./Simulation/ER_Effects_OSEM.RData")
load("./Simulation/ER_Effects_param.RData")

# Check computations
# for (i in 1:500){
# if (!(identical(round(unlist(Effects_OSEM[[i]]),2),round(unlist(Effects_OSEM_dis[[i]]),2)))){
# print(i)
# }
# }

# i = 22,208 not equal, check the position of the mismatch (code below): int 5 on out 8.
# Not relevant, like in pysch one value is 0.0751 the other 0.0750 so the approximation is different
# for (k in 1:20){
# for (j in 1:20){
# if (!(identical(round(unlist(Effects_OSEM[[i]][[k]][[j]]),2),round(unlist(Effects_OSEM_dis[[i]][[k]][[j]]),2)))){
# print(paste0(k,j))
# }
# }
# }

#Generate list of effects for each int-out couple
effects4couple_OSEM<-vector(mode="list", 24)
for (i in 1:n){
  effects4couple_OSEM[[i]]<-vector(mode="list",24)
  for (j in 1:n){
    effects4couple_OSEM[[i]][[j]]<-vector(mode="list",500)
    for (k in c(1:500)){
      effects4couple_OSEM[[i]][[j]][[k]]<-Effects_OSEM[[k]][[i]][[j]]
    }
  }
}

effects4couple_param<-vector(mode="list", 24)
for (i in 1:n){
  effects4couple_param[[i]]<-vector(mode="list",24)
  for (j in 1:n){
    effects4couple_param[[i]][[j]]<-vector(mode="list",500)
    for (k in c(1:500)){
      effects4couple_param[[i]][[j]][[k]]<-Effects_param[[k]][[i]][[j]]
    }
  }
}


effects4couple <- effects4couple_param
### Plot Generation, separate for OSEM and Param: one plot for each variable, showing the effect of 
# switching intervention variable to its lowest to its largest level to all levels of outcome variable
for (i in 1:n){
  png(paste0("ER_total_Int_Param",i,".png"), width = 225, height = 465, units='mm', res = 300)
  lim1<-max(unlist(effects4couple_OSEM[[i]]))
  lim2<-max(unlist(effects4couple_param[[i]]))
  lim=round(max(lim1,lim2),1)
  lim<-round(lim1,1)
  page<-PlotTot(int=i, lim=lim, n=n)
  g<-(wrap_plots(page)+plot_layout(ncol=2))
  print(g)
  dev.off()
}

for (i in 17:17){
  png(paste0("ER_total_Int_New2",i,".png"), width = 1280, height = 720, units='mm', res = 300)
  lim1<-max(unlist(effects4couple_OSEM[[i]]))
  lim2<-max(unlist(effects4couple_param[[i]]))
  lim=round(max(lim1,lim2),1)
  lim<-round(lim1,1)
  page<-PlotTotdouble(int=i, lim=lim)
  g<-(wrap_plots(page)+plot_layout(ncol=5))
  print(g)
  dev.off()
}



# # Plot to compare OSEM and Param Approach
# for (i in 17:17){
#   png(paste0("ER_total_Int_pres",i,".png"), width = 1280, height = 720, units='mm', res = 300)
#   lim1<-max(unlist(effects4couple_OSEM[[i]]))
#   lim2<-max(unlist(effects4couple_param[[i]]))
#   lim=round(max(lim1,lim2),1)
#   lim<-round(lim1,1)
#   page<-PlotTotpres(int=i, lim=c(-lim, lim))
#   g<-(wrap_plots(page)+plot_layout(ncol=5))
#   print(g)
#   dev.off()
# }






###################
# PLOT GRAPHS #####
###################
# True Graph
# Plot DAG
ew <- as.character(round(unlist(edgeWeights(trueDAG)),2))
ew <- ew[setdiff(seq(along=ew), removedEdges(trueDAG))]
names(ew) <- edgeNames(trueDAG)
eAttrs <- list()
eAttrs$label <- ew
attrs <- getDefaultAttrs()
attrs$edge$fontsize <- 8
#plot(trueDAG, edgeAttrs=eAttrs, attrs = attrs)
plot(trueDAG)

# Plot OSEM GRAPH
OSEMfit_point<-results[[2]]$OSEM
g <- as_graphnel(graph_from_adjacency_matrix(OSEMfit_point$DAG))
cpdag_OSEM <- dag2cpdag(g)
#png(paste0("OSEM_CDDAG_ER_Rand",i,".png"), width = 465, height = 225, units='mm', res = 300)
plot(as(cpdag_OSEM, "graphNEL"))
#dev.off()