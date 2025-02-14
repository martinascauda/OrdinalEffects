library(parallel)

cl <- makeCluster(detectCores()-1)


clusterEvalQ(cl,{
  library(BiDAG)
  load("Demo_Sim/ERexample_data.RData")
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
    Paramfit <- getParam(n, trueAM, x,
                               usrpar = list(penType = "other",
                                             L = 5,
                                             lambda = 0))
    Param<- list()
    Param$cuts<-Paramfit[["param"]][["cuts"]]
    Param$S<- Paramfit[["param"]][["Sigma_hat"]]
    
    #Learn both parameters and structure
    OSEMfit <- ordinalStructEM(n, x,
                               usrpar = list(penType = "other",
                                             L = 5,
                                             lambda = 6))
    OSEM<-list()
    OSEM$cuts<-OSEMfit[["param"]][["cuts"]]
    OSEM$S<- OSEMfit[["param"]][["Sigma_hat"]]
    OSEM$DAG <- as.matrix(OSEMfit$DAG)
    
    #Save results  
    res <- list()
    res$Param <- Param
    res$OSEM <- OSEM
    return(res)
  }
})

load("Boot_Sim/BootstrapER_List.RData")

system.time(results <- parLapply(cl,bootstrap.sample,function(x) sim_once(x)))
stopCluster(cl)
save(results, file = "Boot_Sim/Bootstrap_results_ER.RData")