library(parallel)

cl <- makeCluster(detectCores()-1)


clusterEvalQ(cl,{
  library(BiDAG)
  load("Simulation/ER_exampleRand.data.RData")
  source("./OSEMSource/R/OrdinalScore.R")
  source("OrdinalEffects.R")
  insertSource("./OSEMSource/R/spacefns.R",package = "BiDAG")
  insertSource("./OSEMSource/R/usrscorefns.R",package = "BiDAG")
  insertSource("./OSEMSource/R/initpar.R",package = "BiDAG")
  insertSource("./OSEMSource/R/scoreagainstdag.R",package = "BiDAG")
})

clusterEvalQ(cl, {
  sim_once <- function(x) {

    # Learn only parameters (structure known)
    Param<- list()
    prm<-getParam(x,AM)
    Param$cuts<-prm$cuts
    Sigma_hat <- prm$Sigma_hat
    Chol <- getCov(Sigma_hat,AM)
    Param$B <- Chol[[1]]
    Param$Vchol <- Chol[[2]]
    
    
    #Learn both parameters and structure
    OSEMfit <- ordinalStructEM(n, x,
                               usrpar = list(penType = "other",
                                             L = 5,
                                             lambda = 6))
    OSEM<-list()
    OSEM$cuts<-OSEMfit[["param"]][["cuts"]]
    OSEM$S<- OSEMfit[["param"]][["Sigma_hat"]]
    OSEM$DAG <- as.matrix(OSEMfit$DAG)
    res <- list()
    res$Param <- Param
    res$OSEM <- OSEM
    return(res)
  }
})

load("Simulation/BootER_Rand_List.RData")

system.time(results <- parLapply(cl,boot.sample,function(x) sim_once(x)))
stopCluster(cl)
save(results, file = "Bootresults_ER_Rand.RData")